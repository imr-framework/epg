function [Gx,t] = trap2arb(Gxp, system)

total_time = Gxp.riseTime + Gxp.flatTime + Gxp.fallTime;
t = 0:system.gradRasterTime:total_time - system.gradRasterTime;
Npts_flat = ceil(Gxp.flatTime./system.gradRasterTime);

Npts_rise = Gxp.riseTime./system.gradRasterTime;
Npts_fall = Gxp.fallTime./system.gradRasterTime;

Grise =linspace(0,Gxp.amplitude,Npts_rise); 
Gfall = linspace(Gxp.amplitude,0, Npts_fall) ; 
Gflat  = ones(1, Npts_flat).*Gxp.amplitude;

Gx = [Grise Gflat Gfall];
% figure; plot(t,Gx);







