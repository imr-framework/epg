function [ktraj,Grad,trap_len] = gen_kooshball(kmax, N, system)
% Implements the classic Wong and Roos paper, MRM 1994
% Author - Sairam Geethanath

%% Obtain the trajectory first - separated to retain purity of geometry
show =1;
n=1:N;
% alpha = 6; %weighting factor for the poles
% z = atan(alpha*((2*n) - N - 1)/N)./atan(alpha*(N-1)/N);
z = ((2*n) - N - 1)/N; 


phi_t =sqrt(N*pi)*asin(z); %Important to understand this term to avoid SR issues downstream
%% Interpolate phi_t so that diff(phi_t) < 0.05
sin_theta = sqrt(1 - z.^2);

x = cos(phi_t).*sin_theta;
y = sin(phi_t).*sin_theta;

if(length(kmax) >1)
    kzmax = kmax(2);
else
    kzmax = kmax;
end

kx =x*kmax(1); ky =y*kmax(1); kz =z.*kzmax; %Introduced this to make sure kz is also fine
gamma = 42576000;
SRx= get_ktraj2SR(kx, gamma, system);
plot(SRx(2:end));



%% DEBUG only

if(show)
gamma = 42576000;
SRx= get_ktraj2SR([ kx], gamma, system);
SRy= get_ktraj2SR([ ky], gamma, system);
SRz= get_ktraj2SR([kz], gamma, system);

figure(1001); plot(SRx); hold on; plot(system.maxSlew.*ones(size(SRx))./gamma,'k-', 'LineWidth', 2);
hold on; plot(SRy); hold on; plot(SRz); plot(-system.maxSlew.*ones(size(SRx))./gamma,'k-', 'LineWidth', 2);
end

%% Need to interpolate between 0 and the first point for each axis

ktraj = zeros(3,length(kx));
if(show)
figure(1002); plot3(kx,ky,kz,'r*');grid on;
xlabel('kx'); ylabel('ky'); zlabel('kz');
end

ktraj(1,:) = kx; ktraj(2,:) = ky; ktraj(3,:) = kz;
[Gx, Gy, Gz] = k2g(ktraj,system.gradRasterTime);

%% Append gradients to have a modulus of 128 to be compatible with the RF pulses
% SR= get_ktraj2SR(ktraj, gamma, system);
%% Visualize the gradients
dw = system.gradRasterTime;
if(show)
figure(102); 
t = 0:dw:(length(Gx)-1)*(dw);
plot(t,Gx*1e3/gamma);hold on; plot(t,Gy*1e3/gamma); plot(t,Gz*1e3/gamma); grid on; legend('Gx', 'Gy', 'Gz');
xlabel('Time (s)'); ylabel('G (mT/m)');

figure(103);
plot(diff(Gx)/system.gradRasterTime/gamma);hold on; plot(diff(Gy)/system.gradRasterTime/gamma); plot(diff(Gz)/system.gradRasterTime/gamma); grid on; legend('Gx', 'Gy', 'Gz');
xlabel('Time (s)'); ylabel('Slew (T/m/s)');
end
% Npts_proj_grad = Nsp*8; %60us for the pulse + 2 ADC points

%% DEBUG only -  Check the back and forth between ktraj and grad - takes a long time
Gxp.waveform = Gx; Gyp.waveform = Gy; Gzp.waveform = Gz; Gxp.t = [0 dw];
 [kx,ky,kz] = get_g2k(Gxp,Gyp,Gzp);
 if(show)
     figure(1002);hold on;
plot3(kx,ky,kz,'bo');
 end

%% 
Grad.Gx = Gx;
Grad.Gy = Gy;
Grad.Gz = Gz;