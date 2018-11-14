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
% Modified by: Sairam Geethanath, MIRC & CMRRP, fixing SR violations on the
% previous implementations by reinforcing the paper constraints
 
function [kn,Gn,lambda]= vds2D_pulseq_SG(FOV,N,Nshots,alpha,lims)
gamma = 42576000; %Determined from Pulseq - do not change
Gmax=lims.maxGrad/gamma;             %T/m
SRmax=lims.maxSlew/gamma;             %T/m/s
res=FOV/N;               %m
lambda = .5/res; % in m^(-1)
n = (1/(1-(1-Nshots/FOV/lambda)^(1/alpha)));
w = 2*pi*n;

%% Build terms required for the expressions

SRgamma = lims.maxSlew;%Hz/m/s
Ggamma = lims.maxGrad; %Hz/m
Lw = lambda*w;
alpha_p = 0.5*alpha + 1;

%% Generating first interleave
Tes = 1/ ((0.5*alpha + 1)*(sqrt(SRgamma/(Lw*w))));
Tea =1/ ((Ggamma/Lw)*(alpha + 1));

% Ts2a = (Ggamma/((Lw/alpha_p)*Tes^((alpha+1)/alpha_p)))^((alpha +2)/alpha);

% Tea = lambda*w/gamma/Gmax/(alpha+1); % in s
% Tes = sqrt(lambda*w^2/(SRmax*gamma))/(alpha/2+1); % in s
% Ts2a = (Ggamma/((Lw/alpha_p)*Tes^((alpha+1)/alpha_p)))^((alpha +2)/alpha);

Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s
 
if Ts2a<Tes
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+ ...
        ((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a); %This is different
    Tend = Tea;
else
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
    Tend = Tes;
end
 
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));

dt = Tea*1e-4; % in s -- little tricky to figure out
Dt = dt/FOV/abs(k(Tea)-k(Tea-dt)); % in s

t = 0:Dt:Tend; % in s
disp(Dt);
kto = k(t); % in m^-1

%% Derive gradient waveforms 

G = diff(kto)/Dt/gamma; %T/m
Gx = real(G);
Gy = imag(G);

SRx = diff(Gx)/Dt;
SRy = diff(Gy)/Dt;
disp([max(Gx) max(Gy)]);
disp([max(abs(SRx)) max(abs(SRy))]);

[~, dt]= checkgradconst(kto,Dt,SRmax);

%%
t = 0:dt:Tend;
ktest = k(t);
dt=lims. gradRasterTime;
tgrad=0:dt:Tend;
 
%     dist = [0,cumsum(abs(kt(2:end)-kt(1:end-1)))];
    kt_x = interp1(t,real(kto),tgrad,'spline'); 
    kt_y = interp1(t,imag(kto),tgrad,'spline'); 
    kt = complex(kt_x,kt_y);

    
    
    

%%


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
                
 
 
 
 
 


