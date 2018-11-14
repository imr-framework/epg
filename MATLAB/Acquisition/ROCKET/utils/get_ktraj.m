function k = get_ktraj(gx,gy, adc,display)
% gamma_bar = 42.57e6;
% gx_use = gx.amplitude/gamma_bar; %mT/m
% gy_use = gy.amplitude/gamma_bar;
t = linspace(0,gx.flatTime, adc.numSamples);

kx_pre = 0.5.*gx.amplitude.*gx.riseTime; %area of the triangle 
kx = gx.amplitude.*t + kx_pre;

ky_pre = 0.5.*gy.amplitude.*gy.riseTime; %area of the triangle 
ky = gy.amplitude.*t + ky_pre;

k = complex(kx,ky);

if(display)
    figure(1001); grid on; plot(kx,ky,'ro-'); hold on;
    xlabel('kx'); ylabel('ky');
end
