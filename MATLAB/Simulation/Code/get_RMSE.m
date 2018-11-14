function [RMSE_map, Qtri_map, Qimp_map,NRMSE_map] = get_RMSE(Qavg,Qavg_df,Mass_cell)
RMSE_map = zeros(Qavg.dimx, Qavg.dimy,Qavg.dimz);
Qtri_map = RMSE_map;
Qimp_map = RMSE_map;
NRMSE_map = RMSE_map;
R = squeeze(Qavg.index(1:3,:));

%% TODO parallelize this loop - observed offset of 1,1,1 in imp matrix
for k =1:length(R);
    if(Mass_cell(R(1,k),R(2,k),R(3,k))>0)
        RMSE_map(R(1,k),R(2,k),R(3,k)) = RMSE(squeeze(Qavg.avg(k,:,:)), squeeze(Qavg_df(R(1,k)+1,R(2,k)+1,R(3,k)+1,:,:)));
        S = squeeze(Qavg.avg(k,:,:));
        T = squeeze(Qavg_df(R(1,k)+1,R(2,k)+1,R(3,k)+1,:,:));
        Qtri_map(R(1,k),R(2,k),R(3,k)) = sum(abs(S(:)))/numel(S);
        Qimp_map(R(1,k),R(2,k),R(3,k)) = sum(abs(T(:)))/numel(T);
        NRMSE_map(R(1,k),R(2,k),R(3,k))= RMSE_map(R(1,k),R(2,k),R(3,k))./(max(abs(S(:))) - min(abs(S(:)))); 
    end
    
end

function E = RMSE(Q1, Q2)
diff = (Q1 - Q2).^2;
E = sqrt(sum(abs(diff(:)))/numel(Q1));