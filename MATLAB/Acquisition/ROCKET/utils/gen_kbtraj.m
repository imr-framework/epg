 function [k] = gen_kbtraj(Nprojs, Npts, kmax, Tprep,system)
% Gmax=0.9.*lims.maxGrad;

theta = linspace(0, 359, Nprojs);
x = kmax.*cosd(theta);
y = kmax.*sind(theta);
z = kmax.*linspace(0,1, Nprojs);

%% Append points to ensure that they start from 0
Nprep = ceil(Tprep./system.gradRasterTime);
x = [linspace(0,x(1),Nprep)     x(2:end)];
y = [zeros(1,Nprep-1) y];
z = [zeros(1,Nprep-1) z];

k.kx = zeros(Nprojs, Npts);
k.ky= k.kx; 
k.kz = k.kx;

for proj = 1: Nprojs
      k.kx(proj,:) = x(proj).*ones(1,Npts);
      k.ky(proj,:) =  y(proj).*ones(1,Npts);
      k.kz(proj,:) =  z(proj).*ones(1,Npts);  %Use Gz wisely :)
end
%%
k.kx = k.kx(:);
k.ky = k.ky(:);
k.kz = k.kz(:);

plot3(k.kx, k.ky, k.kz);












