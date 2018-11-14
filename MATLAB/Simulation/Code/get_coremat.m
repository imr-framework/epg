function [Bstar,ind_sort,vop_ind,myu_def] = get_coremat(Qinds,ind)
%% Implementation of steps 1 and 2 from the VOP paper
% Step 1 - Identify B*
% Step 2 - Sort matrices
myu_per =0.01; %change back to 0.05

%% Size of Q changes on every iteration.
B = zeros(1,length(ind));
Lambdamin = B;
Qind = squeeze(Qinds(ind,:,:));

%% Identify Bstar Cant use pair-wise to use parforloop
parfor k=1:length(B)
    Qtemp = squeeze(Qind(k,:,:));
    B(k) = norm(Qtemp,2);     
end
[maxB,maxk] = max(B(:)); 
vop_ind = ind(maxk);
myu_def = myu_per* maxB;
Bstar = squeeze(Qind((maxk),:,:));

%% Generate diff matrices, try to make this part smarter
parfor k=1:length(B)
    Qtemp = Bstar - squeeze(Qind(k,:,:));
    Lambdamin(k) = min(eig(Qtemp));
end
[~,ind_sort]=sort(Lambdamin,'descend');
ind_sort = ind(ind_sort);

