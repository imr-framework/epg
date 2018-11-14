%% This file describes the starter version of Tailored MRF implementation through Pulseq
%Author - Sairam Geethanath
%Version - 00
%Date of Version - 4th June 2018
%Version history - Tested playing random flip angles and TR with the tested
%16 shot VDS

%% Create a 2D spiral sequence and export for execution
% This file would serve Marina's project for whole body rapid
% imaging/multiple spiral projects
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
%
% # Create slice selective RF pulse for imaging.
% # Create readout gradients in both directions
% # Loop through number of projections
% # Write the sequence to an open file format suitable for execution on a % scanner.
% Author: Sairam Geethanath
%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assignced default values.
% Pulseq_dir = uigetdir('','Pick the sequences directory');
addpath(genpath('.'));
% system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
% 'MaxSlew',130,'SlewUnit','T/m/s');
gamma = 42576000; % in Hz %Determined from Pulseq - do not change
Gmax = 32; %mT/m
SRmax = 130;%T/m/s
system = mr.opts('MaxGrad',Gmax,'GradUnit','mT/m',...
 'MaxSlew',SRmax,'SlewUnit','T/m/s', 'gradRasterTime', 10e-6, 'rfRingdownTime',10e-6);
%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);
%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3; %One parameter to test
Nx=128; Ny=256;
sliceThickness=5e-3;
load ('TR_brain.mat');
TR_all=TR_all.*1e-3;
TR_all = TR_all - mod(TR_all,1e-3);
TE = 5e-3;%minimum
dx = fov/Nx;
dy =dx;
%% trajectory specific code - design/grad waveforms
Nshots = 16;
alpha = 3; %Arch - resolution vs Rapid - VDS CS
% [ktraj, G] = get_spiral_waveforms();
[ktraj,G,lambda]= vds2D_pulseq_v1(fov,Nx,Nshots,alpha,system);
ktraj = ktraj.*1e-3; %for trajectory to fit into GPI +/-0.5
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2);
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);
%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences.
%
% First, a slice selective RF pulse (and corresponding slice gradient) can % be generated using the |makeSincPulse| function.
%

adc = mr.makeAdc(length(G),'Dwell',system.gradRasterTime);

%% Add 16 different RF pulses 
load('FA_16shots.mat');
for itr = 1:length(FA_16shots)
  [rf, gz] = mr.makeSincPulse(deg2rad(FA_16shots(itr)),system,'Duration',1e-3+0.05e-3.*itr,...
 'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);  
    rf_all(itr) = rf;
    gz_all(itr) = gz; 
end

%% Need to determine multi-slice frequency offsets
Nslices=1;
deltaz=Nslices*sliceThickness;
z=-(deltaz/2):sliceThickness:(deltaz/2);
%% Spoilers
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used.
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',1e-3);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*preTime);
gx = mr.makeArbitraryGrad('x', squeeze(real(G(:,1))),system);
%% Define sequence blocks
% Next, the blocks are put together to form the sequence
total_acq = 128;
numfullim = length(TR_all(1:total_acq))/Nshots; %has to be an integer
% flip=15*pi/180;
% [rf, gz] = mr.makeSincPulse(flip,system,'Duration',1.5e-3,...
%  'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
for slice=1:Nslices
 %freqoffset=gz_all.amplitude.*z(slice);
 %rf_all.freqOffset=freqoffset;

 for ds = 1: numfullim
                 for ns=1:Nshots
                             ind = (ds -1)*Nshots + ns;
                            % disp(ind);
                             TR = TR_all(ind);
                             %% Calculate timing
                             delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf_all(ns))/2);
                             delayTE = delayTE - - mod(delayTE, 1e-5);
                             % delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
                             % - mr.calcDuration(gx) - mr.calcDuration(gzSpoil) - delayTE;
                             delayTR=TR - TE - mr.calcDuration(gx);
                             delayTR = delayTR - - mod(delayTR, 1e-5);
                              
                             delay1 = mr.makeDelay(delayTE);
                             delay2 = mr.makeDelay(delayTR);


                             seq.addBlock(rf_all(ns),gz_all(ns));
                             seq.addBlock(gzReph);
                             gx = mr.makeArbitraryGrad('x', squeeze(real(G(:,ns))),system);
                             gy = mr.makeArbitraryGrad('y', squeeze(imag(G(:,ns))),system);
                             seq.addBlock(delay1);
                             seq.addBlock(gx,gy,adc);
                            %  seq.addBlock(gzSpoil); % removing the spoiler
                             seq.addBlock(delay2);
                 end
 end

end
%% Display the sequence
% ktraj = ktraj./max(abs(ktraj(:)))./2;
%figure(1002);
% seq.plot('TimeRange',[2.7 2.7+18*TR]);
fname = [date,'_Spiral_2DMRF_A3varTR',num2str(Nshots),'_',num2str(alpha),'_',num2str(system.maxSlew/gamma),'_',num2str(Gmax),'_',num2str(fov*1e3),'_'];
save(fname, 'ktrajs' );
%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
% fname = ['Spiral_2D_', num2str(Nshots),'_',num2str(TE), '.seq'];
seq.write([fname,'.seq']);
scan_time = TR.*Nslices.*Nshots;
% disp(scan_time);