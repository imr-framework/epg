

%% load VOPS and local Q matrices
[Filename,Pathname]=uigetfile('*.mat','Pick the VOP matrices');
load(fullfile(Pathname,Filename));


[Filename,Pathname]=uigetfile('*.mat','Pick the local Q matrices');
load(fullfile(Pathname,Filename));

model ='Implemented';

switch model
    case 'Implemented'
%% Reduce dimensions for localQ matrices
[M,N,P,~,~]=size(Qavg.imp);
Qavg.imp = reshape(Qavg.imp,[M*N*P,8,8]);
S = abs(Qavg.imp)>0;
ind = find(squeeze(S(:,4,4)));
Qinds = squeeze(Qavg.imp(ind,:,:));
    case 'Philips'
        Qinds = squeeze(Qavg.avg);
end

%%
matlabpool local 8;
NrofCoils = 8; %NrofCoils
rand_limit=500; %Nrof tries
SAR_VOP = zeros(1,rand_limit);
SAR_Q = SAR_VOP;



%% Check SAR for different waveforms, static shimming to make a generic statement about SAR.
disp('Starting SAR calculations....');
for num_shim_tries=1:rand_limit
    I = complex(rand(NrofCoils,1),rand(NrofCoils,1)); %prob of using a coil out of the 8 is uniform.
    I = I./norm(I); %to make sure norm(I) =1d
    
    [SAR_VOP(num_shim_tries)] = abs(get_SAR(VOPm,I));
    disp(num_shim_tries);
     
     tic;
     [SAR_Q(num_shim_tries)] = abs(get_SAR(Qinds,I)); % takes 14 seconds approx. on 8 cores.
     toc;
end

matlabpool close;

%% Display - Statistical correlation - Figure 8 from Eichfelder
figure;plot(SAR_Q,SAR_VOP,'r*');axis([0 1  0 1]);
hold on;
lower =[0 SAR_Q]; %perfect correlation -> m =1;
upper = lower + myu_def; %over-estimation
plot([0 SAR_Q],lower,'k');
plot([0 SAR_Q],upper,'k');
xlabel('Maximum 10g SAR using exhaustive search');
ylabel('Maximum 10g SAR using VOP');




