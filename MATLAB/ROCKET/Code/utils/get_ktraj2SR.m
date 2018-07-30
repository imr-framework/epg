function indk = get_ktraj2SR(k, gamma, system, show)

dt = system.gradRasterTime;
%% kx
indk =0;viol=0;
kindx_SRb=0;  kindy_SRb =0;  kindz_SRb=0;

kx = squeeze(k(1,:));
Gnx = diff(kx)./dt./gamma; %T/m
SRx = diff(Gnx)./dt;  %T/m/s
L = round(length(SRx)/2);
SRxviol = squeeze(SRx(1:L)) > system.maxSlew./gamma;%exploiting symmetry
if(sum(SRxviol ) >0)
    L = round(length(SRx)/2);
    [ind_vio] = find(SRxviol);
    kindx_SRb = ind_vio(end) + 2; %compensating for the double differential
    viol=1;
end

%%
ky = squeeze(k(2,:));
Gny = diff(ky)./dt./gamma; %T/m
SRy = diff(Gny)./dt;  %T/m/s
SRyviol = squeeze(SRy(1:L)) > system.maxSlew./gamma; %exploiting symmetry
if(sum(SRyviol ) >0)
    [ind_vio] = find(SRyviol);
    kindy_SRb = ind_vio(end) + 2; %compensating for the double differential
    viol=1;
end

%%
kz = squeeze(k(3,:));
Gnz = diff(kz)./dt./gamma; %T/m
SRz = diff(Gnz)./dt;  %T/m/s
SRzviol = squeeze(SRz(1:L)) > system.maxSlew./gamma; %exploiting symmetry
if(sum(SRzviol ) >0)
    [ind_vio] = find(SRzviol);
    kindz_SRb = ind_vio(end) + 2; %compensating for the double differential
        %     kindx_SRe = length(kx) - kindx_SRb;
viol=1;
end


%% Fix these ktrajs
if(viol==1)
[indk] = max([kindx_SRb  kindy_SRb   kindz_SRb]);
end






%% Visualize
if(show)
        figure(112); 
        plot(Gnx*1e3); hold on; plot(Gny*1e3); plot(Gnz*1e3);
        plot(system.maxGrad.*ones(1,length(SRx)).*1e3./gamma, 'r','LineWidth',3);
        plot(-system.maxGrad.*ones(1,length(SRx)).*1e3./gamma, 'r','LineWidth',3);
        ylabel('G(mT/m)'); grid on;

        figure(113); 
        plot(SRx); hold on; plot(SRy); plot(SRz);
        plot(system.maxSlew.*ones(1,length(SRx))./gamma, 'r','LineWidth',3);
        plot(-system.maxSlew.*ones(1,length(SRx))./gamma, 'r','LineWidth',3);
        ylabel('(T/m/s)'); grid on;
end
