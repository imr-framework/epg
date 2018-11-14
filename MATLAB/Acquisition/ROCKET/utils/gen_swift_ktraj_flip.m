function [k,Nprojh] = gen_swift_ktraj_flip(Np, Nproj)

%% Lets fix z and y, determine x - Needs to be a FID trajectory
r = linspace(0,0.5,Np); %Number of points on the radial spoke - 1/2 on each side
Nth = round(sqrt(Nproj));
Nph = round(Nproj./Nth);
Nproj = Nth*Nph;
%%
theta = linspace(-pi/2,pi/2,Nth);
phi = linspace(-pi,pi,Nph);
ind =1;
Nprojh = Nth.*Nph;
kx = zeros(1,Nprojh);
ky =kx; kz=kx;

    for th =1:length(theta)
        for ph = 1:length(phi)
            for rind =1:length(r)
            [kx(ind), ky(ind), kz(ind)] = sph2cart(phi(ph), theta(th), r(rind));
            ind = ind +1;
            end
%             plot3(kx(end),ky(end),kz(end), 'r*'); hold on; grid on; pause(0.1);
        end
    end

kx1 = reshape(kx, [Np, Nproj]);
ky1 = reshape(ky, [Np, Nproj]);
kz1 = reshape(kz, [Np, Nproj]);
    
kx1(:, 2:2:end) = flip(squeeze(kx1(:,2:2:end))); %To avoid going back to 0 and have SR issues
ky1(:, 2:2:end) = flip(squeeze(ky1(:,2:2:end)));
kz1(:, 2:2:end) = flip(squeeze(kz1(:,2:2:end)));

kx = kx1(:).';
ky = ky1(:).';
kz = kz1(:).';

k = cat(1,kx,ky,kz);
%% Display for validation only
%  plot3(kx,ky,kz,'ro');grid on; hold on;drawnow;
% figure(2); 
% subplot(311); plot(x); 
% subplot(312); plot(y);
% subplot(313); plot(z);