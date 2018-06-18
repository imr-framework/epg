%% Create a 3D radial FID sequence and export for execution
% Author: Sairam Geethanath - 
% Descritption - Allows the user to program and run a 3D UTE Radial
% sequence on a Pulseq compatible scanner

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assignced default values.
% Pulseq_dir = uigetdir('','Pick the sequences directory');
addpath(genpath('.'));
dw = 10e-6;
rfdt = 10e-6;
rfrt = 10e-6; %One learning that this has to be there
system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s', 'rfRingdownTime', rfrt, 'rfDeadtime',rfdt);

%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);

%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3;   %One parameter to test
Nx=16; Nz = 1;% Will work with lesser for now - change it back to 64
sliceThickness=30e-3;%was 10e-3 before
TR = 20e-3;
TE = 0.5e-3;%minimum 
dx = fov/Nx;
dz  = sliceThickness/Nz;
%% Define spoke angles - polar (theta); azimuthal (phi)
radp = get_radkparams(dz,dx,fov,'3D');
wr_traj =0;%to write kspace trajectory for reconstruction - 1 yes
theta = linspace(0,179,radp.Ntheta); %Polar
phi = linspace(0,359,radp.Nphi); %azimuthal
%% RF and spoiler
flip=15*pi/180;
st = 30e-3;
% [rf, gz] = mr.makeSincPulse(flip,system,'Duration',1e-3,...
%     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
% [rf,gz]=mr.makeBlockPulse(flip, 'Duration',20e-6,'SliceThickness',st, 'system', system);
[~,gz] = mr.makeBlockPulse(pi/30,'Duration',dw, 'SliceThickness',st ,'system', system);
rf_dur = 1*dw;
[rf] = mr.makeBlockPulse(pi/30,'Duration',rf_dur, 'system', system);
%% Gradients
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 6.4e-3;
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime, 'system', system);
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',0.5e-3);
% gxSpoil = mr.makeTrapezoid('x',system,'Area',gz.area*2,'Duration',3*preTime);
% gySpoil = mr.makeTrapezoid('y',system,'Area',gz.area*2,'Duration',3*preTime);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',7*preTime); %changed this from 7*preTime

%% Calculate timing
% delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf)/2);
delayTE=TE - (mr.calcDuration(rf)/2);
delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ... %Tricky here, not consistent
    - mr.calcDuration(gx) - mr.calcDuration(gzSpoil);

delayTRps = 0.8.*delayTR;%to avoid gradient heating
delayTRps = delayTRps - mod(delayTRps, 1e-5);
delayTRs = 0.2.*delayTR;
delayTRs = delayTRs - mod(delayTRs, 1e-5);


delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTRps);
delay3 = mr.makeDelay(delayTRs);


%% Define sequence blocks
ktraj = zeros(radp.Ns, adc.numSamples,3);
np = 0;

%%
for nt=1:radp.Ntheta
     for ph = 1: radp.Nphi
            % Excitation

             seq.addBlock(rf);
            % Estimate required kspace extents
            kWidth_projx = kWidth.*sind(theta(nt)).*cosd(phi(ph));
            kWidth_projy = kWidth.*sind(theta(nt)).*sind(phi(ph));
            kWidth_projz = kWidth.*cosd(theta(nt));

            % Determine gradient waveforms for each direction
            gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth_projx,'FlatTime',readoutTime);
            gy = mr.makeTrapezoid('y',system,'FlatArea',kWidth_projy,'FlatTime',readoutTime);
            gz = mr.makeTrapezoid('z',system,'FlatArea',kWidth_projz,'FlatTime',readoutTime);
            
            % Determine the kspace trajectory and store it for
            % reconstruction
            G.gx = gx;G.gy = gy;G.gz = gz;
            np = np+1;ktraj(np,:,:) = get_ktraj_v1(G,adc,0);
           
%             seq.addBlock(delay1);
            seq.addBlock(gx,gy,gz,adc);
            seq.addBlock(delay2);
%             seq.addBlock(gzSpoil);%             seq.addBlock(gxSpoil, gySpoil, gzSpoil);
%             delay3 = mr.makeDelay(delayTRs);
            
     end
      disp(nt/radp.Ntheta);
end
%% Normalize the trajectory  and display the sequence
ktrajs = ktraj./max(abs(ktraj(:)))./2;
seq.plot_sg('TimeRange',[0 TR]);


%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
fname = ['Rad3D_FID_v1_nospoil_ noTEdel_noTRdel', num2str(radp.Ns),'_',num2str(Nx),'_',num2str(rf_dur) '.seq'];
seq.write(fname)
if(wr_traj)
     uisave('ktrajs', 'Ktraj' );
end
%%
% % Display the first few lines of the output file
% s=fileread('Rad2D_FID.seq');
% disp(s(1:309))
