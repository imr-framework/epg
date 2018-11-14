%% Comparing inhouse values to procedure by Philips for B1 simulations

% Inhouse B1 fields



% Philips
dirname =uigetdir();

%%
coilb1 = zeros(120,90,386,8);% Hardcoding is ok because this is one model specific.
for k=1:8
coilb1(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_bplus.df']));
end

%%
Bx= df_read(fullfile(dirname,'multix_coil_bxiso.df'));
By= df_read(fullfile(dirname,'multix_coil_byiso.df'));

Bplus = 0.5.*(Bx + 1i*By);
figure;imagesc(abs(Bplus(:,:)));


%% Determine normalization factor per channel.
coilb1_bplus_pre = abs(squeeze(coilb1(:,:,1)));

%% 
 center_field = squeeze(abs(Bplus(:,:,100)));
norm_factor = 1e-6/mean(abs(center_field(:)));






% %%
% ph = pi/2;
% tx_ph = exp(1i*ph);
% Bplus_new = Bplus.*tx_ph;
% figure;imagesc(abs(Bplus_new(:,:)));
% 
% %%
% Bplus2 = circshift(Bplus, [0 0 2]);
% figure;imagesc(abs(Bplus2(:,:)));