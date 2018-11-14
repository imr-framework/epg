%% load VOPS and local Q matrices
[Filename,Pathname]=uigetfile('*.mat','Pick the VOP matrices');
load(fullfile(Pathname,Filename));


[Filename,Pathname]=uigetfile('*.mat','Pick the local Q matrices');
load(fullfile(Pathname,Filename));

%% Reduce dimensions for localQ matrices
[M,N,P,~,~]=size(Qavg.imp);
Qavg.imp = reshape(Qavg.imp,[M*N*P,8,8]);
S = abs(Qavg.imp)>0;
ind = find(squeeze(S(:,4,4)));
Qinds = squeeze(Qavg.imp(ind,:,:));


%%
matlabpool local 8;
NrofCoils = 8; %NrofCoils
x=ones(8,1);

ph = 0:-pi/4:-7*(pi/4);
tx_ph = exp(1i*ph).';
x = x.*tx_ph;

 SAR_Q = calc_SAR(Qinds,x);
 SAR_VOP = calc_SAR(VOPm,x);

SAR_Qscanner = abs(get_SAR(Qinds,x));
SAR_VOP_scanner = abs(get_SAR(VOPm,x));
%% Display - Statistical correlation - Figure 8 from Eichfelder
figure;plot(SAR_Q,SAR_VOP,'r*');axis([0 100  0 100]);
hold on;
lower =[0 SAR_Q]; %perfect correlation -> m =1;
upper = lower + myu_def; %over-estimation
plot([0 SAR_Q],lower,'k');
plot([0 SAR_Q],upper,'k');
xlabel('Maximum 10g SAR using exhaustive search');
ylabel('Maximum 10g SAR using VOP');
matlabpool close;