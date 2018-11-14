function [kx_app, ky_app, kz_app] = append_ktraj(kx, ky, kz, system)

if(nargin<4)
    system.maxSlew = 100; %T/m/s out of 130
   system.maxGrad = 0.30; %T/m out of 0.32
    system.gradRasterTime = 10e-6;% The units is not going to hurt here as the Gmax/SRmax cancel out the factors
end

dt = system.maxGrad/system.maxSlew;

k(1,:) = [0 kx]; k(2,:) = [0 ky]; k(3,:) = [0 kz];
[gx, gy, gz] = k2g(k,dt);
intfact = round(dt/system.gradRasterTime);
Gxp = interp(gx,intfact);
Gyp = interp(gy,intfact);
Gzp = interp(gz,intfact);


[kx_app,ky_app,kz_app] = get_g2k(Gxp,Gyp,Gzp);







