%% Spiral trajectory
% Author : Pavan Poojar 
%% variable declarations
FOV=256e-3;    % FOV in x direction (m)
alpha=3; %density factor (1- Archimedean; >1 - VDS)
Nshots=48;
N=128;         % Number of points in x direction
%% Function to get spiral trajectory - 1 shot
 
[k,tau]= vds2D(FOV,N,Nshots,alpha);
%% Multishot spiral
for s=1:Nshots
    kn(:,s)=k*exp(2*pi*1i*s/Nshots);
    
    plot(real(kn(:,s)),imag(kn(:,s)));hold on;
end
 
 


