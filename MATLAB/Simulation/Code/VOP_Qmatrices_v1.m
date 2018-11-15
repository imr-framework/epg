%% This script describes the implementation of the VOP paper by Eichfelder
clear all;
[Filename,Pathname]=uigetfile('');
 load(fullfile(Pathname,Filename));
 Qavg_df = Qavg.imp;
 matlabpool local 8;
[M,N,P,~,~]=size(Qavg_df);
Qavg_df = reshape(Qavg_df,[M*N*P,8,8]);
clear Qavg;

%% Reduce data size - verify
S = abs(Qavg_df)>0;
ind = find(squeeze(S(:,4,4)));
Qinds = squeeze(Qavg_df(ind,:,:));
Qind = Qinds; %Store Qind for later.
obs_pts = length(ind);
cluster=zeros(1,100); %preallocating for speed.- more than 100 clusters is not effecient
eph = cluster;

%%
vop_map = zeros(M*N*P,1);
indr = 1:length(ind);
figure;
dbstop if error;
VOPm = zeros(500,8,8);

%% Clustering process
VOP=0;label=0;
tic;
while(obs_pts > 0)
    VOP = VOP+1; disp(VOP);
    %% Selection of core matrix
    
    if(VOP==1)
        [Bstar,ind_sorta,myu_def] = get_coremat(Qind,indr);
    else
        [Bstar,ind_sorta] = get_coremat(Qind,indr);
    end
    
    A = Bstar;
    
    if(numel(ind_sorta)>2)
        ind_sort = squeeze(ind_sorta(2:end));
    else
        disp('Reached end of clustering, exiting....');
        break;
    end
    
    %% Set-up problem to find Z*
 
    l=1;e=0;
    Z=zeros(8); %start with a init Z every time for a new cluster.
    cluster_done=0;
    while (cluster_done==0)
       
        % Spectral decomposition 
        Q = A - squeeze(Qind(ind_sort(1),:,:));
        [V,E] = eig(Q);
        Ep = E;
        Ep(Ep<0)=0;
        Em = Ep - E;
        Z_new = V*Em*V';                %Qm = V*E-*V';
        Z = Z + Z_new;
        myu_calc = norm(Z,2);

                 if(myu_calc >= myu_def)
                        if(l>1)
                                 disp([cluster(VOP) myu_calc]);
                                 cluster_done=1;
                                 indr=squeeze(setdiff(indr,(squeeze(ind_sorta(1:l-1)))) );
                                 cluster(VOP)=l-1;


                                 eph(VOP) = -min(eig(A - squeeze(Qind(ind_sorta(l-1),:,:))));
                                 VOPm(VOP,:,:) =A;
                                 label=label+1;
                                 vop_map(ind(ind_sorta(1:l-1)))=label;
                                 S = reshape(vop_map,[M,N,P]);imagesc(abs(squeeze(S(:,45,:))));drawnow;
                        else
                                 disp([cluster(VOP) myu_calc]);
                                 cluster_done=1;
                                 indr=squeeze(setdiff(indr,(squeeze(ind_sorta(1)))) );
                                 cluster(VOP)=l-1;


                                 eph(VOP) = -min(eig(A - squeeze(Qind(ind_sorta(l),:,:))));
                                 VOPm(VOP,:,:) =A;
                                 label=label+1;
                                 vop_map(ind(ind_sorta(1)))=label;
                                 S = reshape(vop_map,[M,N,P]);imagesc(abs(squeeze(S(:,45,:))));drawnow;
                      
                        end
                     
                 else
                    
                    obs_pts = obs_pts -1;
                    lenind = length(ind_sort);
                    switch lenind
                        case 2
                            ind_sort = squeeze(ind_sort(end));
                        case 1
                         disp([cluster(VOP) myu_calc]);
                         cluster_done=1;cluster(VOP)=l;
                         
                         indr=squeeze(setdiff(indr,(squeeze(ind_sorta(1:l)))) );
                         eph(VOP) = -min(eig(A - squeeze(Qind(ind_sorta(l),:,:))));
                         
                         VOPm(VOP,:,:) =A;
                         label=label+1;
                         vop_map(ind(ind_sorta(1:l)))=label;
                         S = reshape(vop_map,[M,N,P]);imagesc(abs(squeeze(S(:,45,:))));drawnow;
                             break;
                        otherwise
                            ind_sort = squeeze(ind_sort(2:end));
                    end
                    A = Bstar + Z;
                    l=l+1;
                    if(mod(l,1000)==0)
                        disp([(l/1e4) myu_calc]);
                    end
%                   
                 end

    end
     
   
    
    
end
 VOPm = squeeze(VOPm(1:VOP-1,:,:));
toc