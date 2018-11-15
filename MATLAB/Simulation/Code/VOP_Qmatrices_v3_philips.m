%% This script describes the implementation of the VOP paper by Eichfelder based on the calculations of Q matrices -imp
clear all;
[Filename,Pathname]=uigetfile('');
load(fullfile(Pathname,Filename));
 Qind = Qavg.avg;
 matlabpool local 8;
clear Qavg Filename Pathname;

%% Reduce data size - verify
obs_pts = size(Qind,1);
cluster=zeros(1,100); %preallocating for speed.- more than 100 clusters is not effecient
% eph = cluster;
normplot = cluster;
vop_ind = cluster;
%%
indr = 1:obs_pts;
ind = indr;
figure;
dbstop if error;
VOPm = zeros(500,8,8);
VOP=0;

tic;
while(obs_pts ~= 0)
    %% Selection of core matrix
    lenindr = length(indr);
    switch lenindr
        case length(ind)
            [Bstar,ind_sorta,vopin,myu_def] = get_coremat(Qind,indr);
%             myu_def=4.9151; %to compare with VOP_imp
        case 1
            disp('Done clustering for all obs pts');break;
            %TODO - if last point is a new cluster
        otherwise
        [Bstar,ind_sorta,vopin] = get_coremat(Qind,indr);
    end
    q=2;
    A = Bstar;
       %% Set-up problem to find Z*
    Z=zeros(8); %start with a init Z every time for a new cluster.
    cluster_done=0;
    obs_pts = obs_pts -1; %corresponding to A.
    while (cluster_done==0)
       
        % Spectral decomposition 
        Q = A - squeeze(Qind(ind_sorta(q),:,:));
        [V,E] = eig(Q);
        Ep = E;
        Ep(Ep<0)=0;
        Em = Ep - E;
        Z_new = V*Em*V';                %Qm = V*E-*V';
        Z = Z + Z_new;
        myu_calc = norm(Z,2);

        
                 if(myu_calc >=myu_def)
                     %% end of current cluster
                     cluster_done=1;
                     VOP = VOP +1;
                     VOPm(VOP,:,:) =A;
                     cluster(VOP) = q-1;%this one broke it, so count upto previous one.
%                      eph(VOP) = -min(eig(A - squeeze(Qind(ind_sorta(q-1),:,:))));
                     normplot(VOP) = norm(Bstar);
                     indr = squeeze(setdiff(indr,(squeeze(ind_sorta(1:q-1))))); 
                     obs_pts_check = sum(cluster) + length(indr);
                     disp([VOP myu_calc obs_pts_check/1e4]);
                     
                     vop_ind(VOP) = ind(vopin);
%                      vop_map(ind(ind_sorta(1:q-1)))=normplot(VOP);
%                      S = reshape(vop_map,[M,N,P]);imagesc(abs(squeeze(S(:,45,:))));drawnow;
                 else
                     %% continue clustering
                     if(q < length(ind_sorta))
                        A = Bstar + Z;
                        obs_pts = obs_pts -1;
                        q = q+1; %for next Q matrix
                                    if(mod(q,1e4)==0)
                                        disp(q/1e4);
                                    end
                      elseif(q==length(ind_sorta))
                         disp('Reached end of clustering process');
                         obs_pts = obs_pts -1;
                         cluster_done=1;
                         VOP = VOP +1;
                         VOPm(VOP,:,:) =A;
                         cluster(VOP) = q;%this is the last one, so here it ends.
%                          eph(VOP) = -min(eig(A - squeeze(Qind(ind_sorta(q),:,:))));
                          normplot(VOP) = norm(Bstar);
                          vop_ind(VOP) = ind(vopin);
%                           vop_map(ind(ind_sorta(1:q-1)))=normplot(VOP);
                         indr = squeeze(setdiff(indr,(squeeze(ind_sorta(1:q))))); 
                         obs_pts_check = sum(cluster) + length(indr);
                         disp([VOP myu_calc obs_pts_check/1e4]);
                         indr=0;
                         break;
                     end
                 end
    end
end
toc;
VOPm = squeeze(VOPm(1:VOP,:,:));
normplot = squeeze(normplot(1:VOP));
vop_ind = squeeze(vop_ind(1:VOP));
matlabpool close;
%% Prepare VOP for writing to scanner compatible format
% VOP_imp = zeros(M,N,P,8,8);
% for k=1:VOP
% [x,y,z]    = ind2sub([M,N,P],vop_ind(k));
% VOP_imp(x,y,z,:,:) =     squeeze(VOPm(k,:,:));
% end

%% Plot norm of the VOPs - Figure 5 from Eichfelder
figure;plot(1:size(VOPm,1),abs(normplot),'ko');
xlabel('Index of VOP');ylabel('Spectral norm of VOP');
