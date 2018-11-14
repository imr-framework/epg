function k = get_ktraj_v1(G, adc,display)
% gamma_bar = 42.57e6;
% gx_use = gx.amplitude/gamma_bar; %mT/m
% gy_use = gy.amplitude/gamma_bar;

gx = G.gx;
gy = G.gy;
gz = G.gz;
t = linspace(0,gx.flatTime, adc.numSamples);

kx_pre = 0.5.*gx.amplitude.*gx.riseTime; %area of the triangle 
kx = gx.amplitude.*t + kx_pre;

ky_pre = 0.5.*gy.amplitude.*gy.riseTime; %area of the triangle 
ky = gy.amplitude.*t + ky_pre;

kz_pre = 0.5.*gz.amplitude.*gz.riseTime; %area of the triangle 
kz = gz.amplitude.*t + kz_pre;

k(:,1) = kx; k(:,2) = ky; k(:,3) =kz;


if(display)
    figure(1001); grid on; plot3(kx,ky,kz,'ro-'); hold on;
%    drawnow;
end
 xlabel('kx'); ylabel('ky'); zlabel('kz');