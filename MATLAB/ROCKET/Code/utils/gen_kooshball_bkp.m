function [ktraj,Grad,trap_len] = gen_kooshball(kmax, N, system)
% Implements the classic Wong and Roos paper, MRM 1994
% Author - Sairam Geethanath

%% Obtain the trajectory first - separated to retain purity of geometry
gamma = 42576000;

n=1:N;
zp = (((2*n) - N - 1)/N);
ramp = round(0.1*N);
mid = round(0.8.*N);

zru = linspace(-0.999, -0.7500, ramp);
zrd = linspace(0.7500,0.999,ramp);
zmid = linspace(-0.7499,0.7499,mid);

z = [zru zmid zrd];
figure(101); plot(zp); grid on; hold on; plot(z);
%%
phi_t =sqrt(N*pi)*asin(z); %Important to understand this term to avoid SR issues downstream
plot(phi_t);

phi_t =sqrt(N*pi)*asin(zp); %Important to understand this term to avoid SR issues downstream
plot(phi_t);hold on;
%% SR challenges
SRphi = diff(phi_t);

SRphi_norm = SRphi./max(SRphi(:));
figure; plot(phi_t); hold on; plot(SRphi_norm);
sin_theta = sqrt(1 - z.^2);

x = cos(phi_t).*sin_theta;
y = sin(phi_t).*sin_theta;

if(length(kmax) >1)
    kzmax = kmax(2);
else
    kzmax = kmax;
end

kx =x*kmax(1); ky =y*kmax(1); kz =z.*kzmax; %Introduced this to make sure kz is also fine


SRx= get_ktraj2SR(kx, gamma, system);
SRy= get_ktraj2SR(ky, gamma, system);
SRz= get_ktraj2SR(kz, gamma, system);
show =1;
if(show)
figure; plot(SRx); hold on; plot(system.maxSlew.*ones(size(SRx))./gamma,'k-', 'LineWidth', 2);
hold on; plot(SRy); hold on; plot(SRz); 
end
SR  = max([SRx SRy SRz]);
if(SR > system.maxSlew)
    error('Increase number of points in the trajectory'); %This has to be fixed
end
%% Need to interpolate between 0 and the first point for each axis
% [kx_app, ky_app, kz_app] = append_ktraj(kx(1), ky(1), kz(1));
d = system.prepgradtime; % This is a problem to have this hardcoding here
trap_len = d/system.gradRasterTime;

kxp = linspace(0,kx(1),trap_len);
kyp = linspace(0,ky(1),trap_len);
kzp = linspace(0,kz(1),trap_len);

kx = [kxp kx]; ky = [kyp ky]; kz = [kzp kz];
ktraj = zeros(3,length(kx));
% figure(1002); plot3(kx,ky,kz);grid on;
% xlabel('kx'); ylabel('ky'); zlabel('kz');
ktraj(1,:) = kx; ktraj(2,:) = ky; ktraj(3,:) = kz;
[Grad.Gx, Grad.Gy, Grad.Gz] = k2g(ktraj,system.gradRasterTime);
% SR= get_ktraj2SR(ktraj, gamma, system);
%% Visualize the gradients
dw = system.gradRasterTime;
figure(102); 
t = 0:dw:(length(Gx)-1)*(dw);
plot(t,Gx*1e3/gamma);hold on; plot(t,Gy*1e3/gamma); plot(t,Gz*1e3/gamma); grid on; legend('Gx', 'Gy', 'Gz');
xlabel('Time (s)'); ylabel('G (mT/m)');

figure(105);
plot(diff(Gx)/system.gradRasterTime/gamma);hold on; plot(diff(Gy)/system.gradRasterTime/gamma); plot(diff(Gz)/system.gradRasterTime/gamma); grid on; legend('Gx', 'Gy', 'Gz');
xlabel('Time (s)'); ylabel('Slew (T/m/s)');
% Npts_proj_grad = Nsp*8; %60us for the pulse + 2 ADC points

%% DEBUG only -  Check the back and forth between ktraj and grad - takes a long time
Gxp.waveform = Gx; Gyp.waveform = Gy; Gzp.waveform = Gz; Gxp.t = [0 dw];
 [kx,ky,kz] = get_g2k(Gxp,Gyp,Gzp);
plot3(kx,ky,kz,'bo');