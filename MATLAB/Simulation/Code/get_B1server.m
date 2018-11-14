function B1server = get_B1server(A)
% function B1server = get_B1server(A,Bx,By, Air);
% Calculate B1Server values based on the iso-slice
Ncoils =8;
%%  Read input files - need to make these inputs later
Bx = df_read('multix_coil_bxiso.df');
By = df_read('multix_coil_byiso.df');
Air = df_read('multix_coil_bair_seg.df');
Air(Air ==-1) =0;

ph = 0:pi/4:7*(pi/4);
tx_ph = exp(1i*ph).';
% A = 2.*ones(8,1);

%% Use airseg to exclude voxels laterd
for nc =1:Ncoils
    Bx(:,:,nc) = Bx(:,:,nc).*Air.*tx_ph(nc).*A(nc);
    By(:,:,nc) = By(:,:,nc).*Air.*tx_ph(nc).*A(nc);
end

Bx_sum = sum(Bx,3);
By_sum = sum(By,3);


%% Obtain the mean B1 value experienced due to the shims applied.
B = Bx_sum + 1i*By_sum;
B1server =0.5.* mean(abs(B(B~=0))); % Tesla






