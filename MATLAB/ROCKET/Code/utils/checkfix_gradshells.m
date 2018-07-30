function   [Gx, Gy, Gz] = checkfix_gradshells(Gx,Gy,Gz,system)
%% Check Gmax and SR violations
gamma = 42576000;
Gx_viol = Gx > system.maxGrad;
Gy_viol = Gy > system.maxGrad;
Gz_viol = Gz > system.maxGrad;

SRx_viol = (diff(Gx)./system.gradRasterTime) > system.maxSlew;
SRy_viol = (diff(Gx)./system.gradRasterTime) > system.maxSlew;
SRz_viol = (diff(Gx)./system.gradRasterTime) > system.maxSlew;

%% Visualize the violations if any
figure(112);
plot(Gx/gamma); hold on;
plot(Gy/gamma); hold on;
plot(Gz/gamma); hold on;
plot(system.maxGrad.*ones(1,length(Gx_viol))./gamma,'k', 'LineWidth',2);

%%
figure(113);
plot(diff(Gx)./system.gradRasterTime./gamma); hold on;
plot (diff(Gx)./system.gradRasterTime./gamma) ; hold on;
plot (diff(Gx)./system.gradRasterTime./gamma) ; hold on;
plot(system.maxSlew.*ones(1,length(Gx_viol))./gamma,'k', 'LineWidth',2);