%% Understanding Philips B1 normalization



%% Check definition of Bplus - cropped
% Philips
dirname =uigetdir();


coilb1 = zeros(120,90,386,8);
parfor k=1:8
coilb1_plus(:,:,:,k)= df_read(fullfile(dirname,['multix_coil',num2str(k),'_bplus.df']));
% coilb1(:,:,k) = squeeze(s(:,:,193));
end

Bx= df_read(fullfile(dirname,'multix_coil2_bx.df'));
By= df_read(fullfile(dirname,'multix_coil2_by.df'));
Bplus = 0.5.*(Bx + 1i*By);
Bplus_coil1 = squeeze(coilb1_plus(:,:,:,2));

% 
Diff_definition = Bplus - Bplus_coil1;
disp(sum(abs(Diff_definition(:))));


% Visual comparison
figure;subplot(121);imagesc(abs(Bplus(:,:,193)));
subplot(122);imagesc(abs(coilb1_plus(:,:,193,2)));

 
% DEFINITION CONSISTENT, PROCEED FURTHER!

%% Compare B1plus - with and without norm
dirname =uigetdir();
Bx_norm= df_read(fullfile(dirname,'multix_coil2_bx.df'));
By_norm= df_read(fullfile(dirname,'multix_coil2_by.df'));
Bplus_norm = 0.5.*(Bx_norm + 1i*By_norm);
Bplus_norm_file= df_read(fullfile(dirname,'multix_coil2_bplus.df'));

% Visual comparison
figure;subplot(121);imagesc(abs(Bplus(:,:,193)));
subplot(122);imagesc(abs(Bplus_norm(:,:,193)));
Norm_factor = Bplus_norm_file./Bplus; % 213, 228

% Verify
Bplus_iso = squeeze(Bplus(:,:,193));
disp(abs(mean(abs(Bplus_iso(:)))*Norm_factor(1)));




