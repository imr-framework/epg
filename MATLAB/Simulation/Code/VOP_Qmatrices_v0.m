%% This script describes the implementation of the VOP paper by Eichfelder
% % matlabpool local 8;
% [M,N,P,~,~]=size(Qavg_df);
% Qavg_df = reshape(Qavg_df,[M*N*P,8,8]);

%% Reduce data size - verify
S = Qavg_df>0;
ind = find(squeeze(S(:,4,4)));
Qind = squeeze(Qavg_df(ind,:,:));
obs_pts = length(ind);
%% Clustering process
VOP=0;
while(obs_pts > 0)
    VOP = VOP+1;
    disp(VOP);
    %% Selection of core matrix
    [Bstar,ind_sorta,myu_def] = get_coremat(Qind,ind);
    A = Bstar;
    ind_sort = squeeze(ind_sorta(2:end));
    t0 = cputime;
    %% Set-up problem to find Z*
    
    l=2;e=0;Z=zeros(8); 
    cluster_done=0;
      while (cluster_done==0)
        Q = A - squeeze(Qind(ind_sort(1),:,:));
        [U,E,V] = svd(Q);
        Ep = E;
        if(sum(Ep(:) < 0))
            disp('One here');
            tdiff = cputime -t0;
            disp(tdiff);
        end
        Ep(Ep <0)=0;
        Em = Ep - E;
        Z = U*Em*V';                %Qm = V*E-*V';
      
        if(sum(abs(Z(:)))>0)
            disp('Found Z')
        end
        
        myu_calc = norm(Z,2);
                 
                             if(myu_calc > myu_def)
                                 cluster_done=1;
                                 ind=setdiff(ind, squeeze(ind_sorta(1:l)));
                             else
                                 
                                 A = A + Z;
                                obs_pts = obs_pts -1;
                                        switch obs_pts
                                            case 0
                                                break;
                                           case 1
                                              ind_sort = (ind_sort);
                                           otherwise
                                            ind_sort = squeeze(ind_sort(2:end));
                                        end
                                l=l+1;
                                if(mod(l,1000)==0)
                                    disp(l);
                                end
                             end

        
        
        
        
        
        
        
        
        
      end
    
    
end


%% Collection of Matrices in a cluster



%% Determination of a VOP for a cluster