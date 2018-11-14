%% GenerateSpiralTraj.m
%
%[k]=GenerateSpiralTraj_v1()
% Generates a spiral k-space trajectory with a method adapted from [1]
% (the paper is a bit buggy).
%
% INPUTS:   * field of view in meters
%           * resolution (Nyquist distance) in meters
%           * undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%           * undersampling factor (normalized units, positive, may be lower than one)
%           * number of interleaves
%           * variable density factor (alpha in [1])
%           * maximum amplitude (in T/m)
%           * maximum slew rate (in T/m/s)
%           * resampling trajectory (true or false)
%           * analysis option (true or false)
%
% [1]   "Simple Analytic Variable Density Spiral Design"
%       Dong-hyun Kim, Elfar Adalsteinsson, and Daniel M. Spielman
%       Magnetic Resonance in Medicine 50:214-219 (2003)
%
% Copyright, Matthieu Guerquin-Kern, 2012
% Modified by : Pavan Poojar, MIRC
% Modified by: Sairam Geethanath, MIRC, CMRRP - Fixed value of Dt to ensure
% SRmax compliance - critical for playing out spirals
 
function [kn,Gn,lambda]= vds2D_pulseq(FOV,N,Nshots,alpha,lims)
gamma =  42576000; % in Hz  %Determined from Pulseq - do not change
SRsafety = 0.7;
Gmax=lims.maxGrad/gamma;             %T/m
SRmax=SRsafety.*lims.maxSlew/gamma;             %T/m/s - safety measures
res=FOV/N;               %m
%% Generating first interleave

lambda = .5/res; % in m^(-1)
n = (1/(1-(1-Nshots/FOV/lambda)^(1/alpha)));
w = 2*pi*n;
Tea = lambda*w/gamma/Gmax/(alpha+1); % in s
Tes = sqrt(lambda*w^2/(SRmax*gamma))/(alpha/2+1); % in s
Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s
 
if Ts2a < Tes 
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a);
    Tend = Tea;
else
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
    Tend = Tes;
end
 
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));
%% Choosing the correct dt to avoid SR violation
% dt = Tea*1e-4; % in s -- little tricky to figure out
dt=lims. gradRasterTime;

t =0:dt:Tend;
kuse = k(t);


 [~, dt2] =checkgradconst(kuse,dt,SRmax);


%%  Should start from here again
kn = zeros(length(kt), Nshots);
Gn = zeros(length(kt), Nshots);
            for s=1:Nshots
                kn(:,s)=kt*exp(2*pi*1i*s/Nshots);%m-1
%                 Gn(:,s) =[0  (diff(squeeze(kn(:,s))).*1e3)/dt]; %mT/m
                km(1,:) = real(squeeze(kn(:,s)));
                km(2,:) = imag(squeeze(kn(:,s)));
                [Gx, Gy] = k2g(km,dt);
                Gx(1) =0;
                Gy(1) =0;
                Gn(:,s) = complex(Gx,Gy);%Hz/m
                plot(real(kn(:,s)),imag(kn(:,s)));hold on;
            end


%         if(( max(abs(Gn(:))) - mod(max(abs(Gn(:))),0.1)) > (Gmax*1e3)) %mT/m
%             error('Design error')
%         end

end
                
 
 
 
 
 


