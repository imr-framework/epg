%% Prepare the GOA gradients for waveform 

%Loads G_all which is a cell
SRmax = 90;%T/m/s
Gmax = 0.03; %T/m 
system.gradRasterTime = 10e-6;
len = 0;
for shot = 1:length(G_all)
    G_shot = G_all{shot};
    Gx = real(G_shot); delGx = abs(Gx(1));
    Gy = imag(G_shot); delGy =abs(Gy(1));
    
    tG = (max([(delGx/SRmax) (delGy/SRmax)]));
    Npts = ceil(tG/system.gradRasterTime);
    
    Gx_new  = [linspace(0,Gx(1),Npts) Gx(2:end)];
    Gy_new  = [linspace(0,Gy(1),Npts) Gy(2:end)];
    
    figure(123);plot(diff(Gx_new)./system.gradRasterTime); hold on;
    
    sd = length(Gx_new);
    len = max([len,sd]);
    
    G_Goa{shot} = complex(Gx_new, Gy_new);
    
    figure(121);plot(Gx_new); hold on;
    figure(121);plot(Gy_new);
    
end
% if(mod(len,2)~=0)
%     len = len +1;
% end
% t =0:system.gradRasterTime: (len-1)*system.gradRasterTime;
% 
% G_Goa_store = zeros(len, length(G_all));
% 
% for shot =1 :length(G_all)
%     Gx = real(G_Goa{shot});
%     Gy = imag(G_Goa{shot});
%     tcurr = 0:system.gradRasterTime: (length(Gy)-1)*system.gradRasterTime;
%     Gx_new = interp1(tcurr, Gx, t,'spline','extrap');
%     Gy_new = interp1(tcurr, Gy, t,'spline','extrap');
%     if(isnan(Gx_new))
%         disp('error')
%     end
%     
%     G_Goa_store(:,shot) = complex(Gx_new,Gy_new);
%     
% end
% 
