%% check B fields
dirname=uigetdir('');
filename='BxIso.df';
bx_iso = df_read(fullfile(dirname,filename));
filename='ByIso.df';
by_iso = df_read(fullfile(dirname,filename));

%% Bx and By fields
Bplus = bx_iso + 1i.*by_iso;
Bplus_quad = sum((Bplus),3);
figure;imagesc(abs(Bplus_quad));
figure;imagesc(rad2deg(angle(Bplus_quad)));
%%
figure;imagesc(abs(Bplus(:,:)));

