
%%
dirname = uigetdir('');
parfor k=1:8 %num_channels
    Ex(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ex.df'])); 
    Ey(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ey.df'])); 
    Ez(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ez.df'])); 

end

Ex = sum(Ex,4);
Ey = sum(Ey,4);
Ez = sum(Ez,4);

E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);

%%
figure;imagesc(abs(squeeze(E(30,:,:))));
figure;imagesc(squeeze(abs(E(:,15,:))));