%% Data driven spirals 
% Author: Pavan Poojar
%% Load the .mat files - image and circular mask 
load 'R_Invitro_N.mat'; % Load image
k1=fftshift(ifft2(I));  % Convert to k-space

load 'cir_mask_256_256.mat'; % Load circular mask
mask=~Mdes;
k1=fftshift(ifft2(I)); % k-space 

[M,N] = size(k1);% Read matrix size of image/k-space

%% Thresholded k-space
sft=k1.*Mdes; % Thresholded k-space
k_mean=abs(mean(sft(:)));   % Mean of the thresholded k-space
k_std=abs(std(sft(:)));     % Standard deviation of the thresholded k-space
k=ones(M,N);
for i=1:M
for j=1:N
    if abs(sft(i,j))<k_mean+k_std/2 % threshold to find alpha (Mean+(SD/2))
        k(i,j)=0;
    end
end
end
%% Function to generate a circular mask

[Mdes,alpha] = circle_mask(M,N,k);
%% Pulseq 
%+++++++++++++++++++++++++++++++++++++++++++++
% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assignced default values.
 Pulseq_dir = uigetdir('','Pick the sequences directory');
addpath(genpath('.'));
addpath(genpath(Pulseq_dir));
% system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
%     'MaxSlew',130,'SlewUnit','T/m/s');
gamma =  42576000; % in Hz  %Determined from Pulseq - do not change
Gmax = 32; %mT/m
SRmax = 130;%T/m/s
system = mr.opts('MaxGrad',Gmax,'GradUnit','mT/m',...
    'MaxSlew',SRmax,'SlewUnit','T/m/s', 'gradRasterTime', 10e-6);


%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);
%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3;   %One parameter to test
Nx=256; Ny=256;
sliceThickness=5e-3;
TR = 50e-3;
TE = 5e-3;%minimum 
dx = fov/Nx;
dy  =dx;
Nshots = 16;
alpha = 2; % 

phiNd=360/Nshots;
X = ['Normal angle','=',num2str(phiNd),' deg.'];
disp(X);
phiN=0;
           
%% trajectory specific code - design/grad waveforms
% [ktraj, G] = get_spiral_waveforms();
[ktraj,G,lambda]= vds2D_pulseq_v1_N_GA_TGA(fov,Nx,Nshots,alpha,system,phiN);

ktraj = ktraj.*1e-3; %for trajectory to fit into GPI +/-0.5
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2); 
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);

Nslices=1;
deltaz=Nslices*sliceThickness;
z=-(deltaz/2):sliceThickness:(deltaz/2);
%% Slice selection
flip=15*pi/180;
[rf, gz] = mr.makeSincPulse(flip,system,'Duration',1.5e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
adc = mr.makeAdc(length(G),system,'Dwell',system.gradRasterTime);

%% Need to determine multi-slice frequency offsets

%% Spoilers
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used.
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',1e-3);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*preTime);
gx = mr.makeArbitraryGrad('x', squeeze(real(G(:,1))),system);
%% Calculate timing
delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf)/2);
% delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
%     - mr.calcDuration(gx) - mr.calcDuration(gzSpoil) - delayTE;
delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
    - mr.calcDuration(gx) - mr.calcDuration(gzSpoil);
delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTR);

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
for slice=1:Nslices
    freqoffset=gz.amplitude*z(slice);
    rf.freqOffset=freqoffset;
for ns=1:Nshots
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
       
             gx = mr.makeArbitraryGrad('x', squeeze(real(G(:,ns))),system);
             gy = mr.makeArbitraryGrad('y', squeeze(imag(G(:,ns))),system);
           
            seq.addBlock(delay1);
            seq.addBlock(gx,gy,adc);
            seq.addBlock(gzSpoil);
            seq.addBlock(delay2);
end
end
%% Display the sequence 
% ktraj = ktraj./max(abs(ktraj(:)))./2;

figure(1002);
seq.plot('TimeRange',[TR 2*TR]);
%seq.plot();
%fname = ['Spiral_2D_',num2str(Nshots),'s','_',num2str(alpha),'a','_',num2str(system.maxSlew/gamma),'_',num2str(Gmax),'_',num2str(phiNd),'deg'];
fname = ['Spiral_2D_',num2str(Nshots),'s','_',num2str(alpha),'a','_',num2str(system.maxSlew/gamma),'_',num2str(Gmax)];

%save(fname, 'ktrajs' );
save([fname,'.mat'], 'ktrajs' );
%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
% fname = ['Spiral_2D_', num2str(Nshots),'_',num2str(TE), '.seq'];
%seq.write([fname,'.seq']);



