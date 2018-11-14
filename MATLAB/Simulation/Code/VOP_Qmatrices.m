%% This script describes the implementation of the VOP paper by Eichfelder
matlabpool local 2;
myu = 0.15; %SAR over-estimation parameter
[M,N,P,~,~]=size(Qavg_df);
Qavg_df = reshape(Qavg_df,[M*N*P,8,8]);

%% Reduce data size - verify
S = Qavg_df>0;
ind = find(squeeze(S(:,4,4)));
Qind = squeeze(Qavg_df(ind,:,:));

%% Selection of core matrix
[Bstar,ind_sort] = get_coremat(Qind,ind);

%% Set-up optimization problem for clustering
obs_pts = length(ind);
l=2;e=0;Z=zeros(8);
Cluster = squeeze(ind_sort(1:2));


%% Clustering process
while(obs_pts > 0)
    
    %
    
    
    
    
    
end


%% Collection of Matrices in a cluster



%% Determination of a VOP for a cluster