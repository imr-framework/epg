%% Obtain an estimate for sigma/rho per channel
dirname = uigetdir();
Nc=8;
%%
temp =importdata(fullfile(dirname,['channel_1_Ex.mat'])); 
Ex = zeros(size(temp,1), size(temp,2), size(temp,3), 8);
Ey = Ex;
Ez = Ex;
clear temp;
parfor k=1:8 %num_channels
    Ex(:,:,:,k) =importdata(fullfile(dirname,['channel_',num2str(k),'_Ex.mat'])); 
    Ey(:,:,:,k) =  importdata(fullfile(dirname,['channel_',num2str(k),'_Ey.mat'])); 
    Ez(:,:,:,k) =  importdata(fullfile(dirname,['channel_',num2str(k),'_Ez.mat'])); 
end

%% Generate Epwr per channel
Epwr = zeros(size(Ex));
SbR_hat = ones(size(Ex));
tic;
z=100;
for k=1:Nc
    for x=1:size(Ex,1)-1
        parfor y =1:size(Ex,2)-1
%             tic;
%            parfor z = 1: size(Ex,3)-1
                X = [x y z];
                Epwr(x,y,z,k) = gen_E12ptQ(squeeze(Ex(:,:,:,k)),squeeze(Ey(:,:,:,k)),squeeze(Ez(:,:,:,k)),X,SbR_hat);
%            end
%             toc;
        end
%         if(mod(x,10)==0)
            disp(x);
%         end
    end
end    
%% Trying to speed it up.
% Epwr = zeros(size(Ex));
% SbR_hat = ones(size(Ex));
% tic;
% h=0;
% Xind = zeros(3840000,1);
% X = Xind;
% Y = Xind;
% Z = Xind;
%     for x=1:size(Ex,1)-1
%         for y =1:size(Ex,2)-1
%             for z = 1: size(Ex,3)-1
%                 h = h+1;
%                 Xind(h) = sub2ind(size(Ex),x, y, z);
%                 X(h) = x;
%                 Y(h) =y;
%                 Z(h) =z;
% 
%            end
%         end
%         disp(x);
%     end
% 
%  %%
% Ef = zeros(8,length(Xind));
% for k=1:Nc
%     Ex_temp = squeeze(Ex(:,:,:,k));
%     Ey_temp = squeeze(Ey(:,:,:,k));
%     Ez_temp = squeeze(Ez(:,:,:,k));
%     
%    parfor h=1:length(Xind)
%          Epwr(h) = gen_E12ptQ(Ex_temp,Ey_temp,Ez_temp,X(h),SbR_hat);
%     end
%     disp(k);
%     Ef(k,:) = Epwr;
% end
% 
% % Re-cast
% E_Nc = zeros(size(Ex));
% for k=1:Nc
%     for h=1:length(Xind)
%         E_Nc(X(h),Y(h),Z(h),k) = Ef(k,h);
%     end
%     disp(Nc);
% end
%    
%%
SAR = zeros(121,161,201,8);
for k=1:Nc
    SAR(:,:,:,k) = importdata(['DHC_SAR-el',num2str(k),'.mat']);
end

%% Divide the SAR matrices by E matrices to get an estimate of SigmabyRho
SigmabyRho_hat = squeeze(SAR(:,:,100,:)./Epwr(:,:,100,:));


