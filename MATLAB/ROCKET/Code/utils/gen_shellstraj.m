function [k,G] = gen_shellstraj(FOV, Np, Ns, system)

%% Calculate relevant inputs
delk = 1/FOV; %Assume isotropic
k0 = (0.5*delk):delk:(Ns-0.5)*delk;
p =0:1:Np-1;

%%
for S = 1:Ns
    Ms = round((8/sqrt(3))*pi*(S^2/Np));
    kzp = zeros(Np, Ms); kyp = zeros(Np, Ms);   kxp = zeros(Np, Ms); 
    if(Ms < 1)
        Ms = 1;
    end
       for ms =1:Ms
            kzp(:,ms) = k0(S)*((2*p) - (Np-1))./(Np-1); 
%              kzp(:,ms) = k0(S).* (alpha*atan((2*p) - (Np-1)./(Np)))./atan(alpha*(Np-1)/Np); 
       
            
            X = sqrt(1 - (kzp(:,ms)./k0(S)).^2);
            theta = (sqrt(Np*pi/Ms)*asin(kzp(:,ms)./k0(S))) + (2*ms*pi/Ms);
            
            kxp(:,ms) = k0(S).*X.*cos(theta);
            kyp(:,ms) = k0(S).*X.*sin(theta);
%             figure(111); plot3(squeeze(kxp(:,ms)),squeeze(kyp(:,ms)),squeeze(kzp(:,ms))); hold all; grid on;
            
            %% Generate gradient waveform
            kuse(1,:) = squeeze(kxp(:,ms));  kuse(2,:) = squeeze(kyp(:,ms)); kuse(3,:) = squeeze(kzp(:,ms)); 
            indk= get_ktraj2SR(kuse, system.gamma, system,0);
            
            if(indk >0)
               kfixed = fix_ktraj(kuse,indk,system,system.gamma);
            else
                kfixed = kuse;
            end
            %% for Debug only
            if(mod(S,10)==0)
              get_ktraj2SR(kfixed, system.gamma, system, 1);
            end
            [Gx, Gy, Gz] = k2g(kfixed,system.gradRasterTime);
            
       
             
            
       end
        disp(S);
        %% Store kspace
         kxp_all{S} = kxp; Gxp_all{S} = Gx;
         kyp_all{S} = kyp;  Gyp_all{S} = Gy;
         kzp_all{S} = kzp; Gzp_all{S} = Gz;
        
        %% 
        
end

k.kx = kxp_all;
k.ky = kyp_all;
k.kz = kzp_all;

G.Gx = Gxp_all;
G.Gy = Gyp_all;
G.Gz = Gzp_all;