function [Mdes, FOVx] = circle_mask(M,N,k)
for i=1:0.5:12
FOVx=i;
M=256;
N=256;
radius = FOVx/10;
centerW = M/2; % to generate circle at the center
centerH = N/2;
[W,H] = meshgrid(1:M,1:N);
mask = (sqrt((W-centerW).^2 + (H-centerH).^2))*1e-2 < radius; % generate cicular mask
%picture = 1;
% Mdes = picture .* mask;
%s=k.*Mdes;
Mdes=mask;
 if(sum(Mdes(:))>sum(k(:)));
     break;
 else
     continue;
 end

end