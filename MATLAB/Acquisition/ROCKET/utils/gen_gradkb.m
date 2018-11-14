 function [Gx,Gy,Gz,Nint] = gen_gradkb(Nprojs, Npts, lims, t)
Gmax=0.9.*lims.maxGrad;

theta = linspace(0, 359, Nprojs);
SRsafety = 0.8;
SRmax=SRsafety*lims.maxSlew;

x = Gmax.*cosd(theta);
y = Gmax.*sind(theta);
z = Gmax.*linspace(0,1, Nprojs);

if(nargin < 4)
t = x(1)/SRmax;
end
dt = lims.gradRasterTime;
Nint = ceil(t/dt);

G.Gxprep = linspace(0,x(1),Nint);
G.Gyprep = zeros(size(G.Gxprep));
G.Gzprep = G.Gyprep;

G.Gx = zeros(Nprojs, Npts);G.Gy= G.Gx; G.Gz = G.Gx;

for proj = 1: Nprojs
      G.Gx(proj,:) = x(proj).*ones(1,Npts);
      G.Gy(proj,:) =  y(proj).*ones(1,Npts);
      G.Gz(proj,:) =  z(proj).*ones(1,Npts);  %Use Gz wisely :)
%       figure(101); plot(squeeze(G.Gx(proj,:))); hold on; plot(squeeze(G.Gy(proj,:))); hold on;  plot(squeeze(G.Gz(proj,:)));
end
G.Gx = G.Gx.'; G.Gy = G.Gy.'; G.Gz = G.Gz.'; 
Gx =[G.Gxprep(:); G.Gx(:)];Gy = [G.Gyprep(:); G.Gy(:)];Gz = [G.Gzprep(:) ;G.Gz(:)];


%% Append gradients to ensure that they start from 0










