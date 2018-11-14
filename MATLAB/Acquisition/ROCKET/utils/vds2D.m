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
 
function [k,tge,lambda]= vds2D(FOV,N,Nshots,alpha)

Gmax=33e-3;%T/m
SRmax=100 ; %T/m/s
res=FOV/N;
%% Generating first interleave
gamma = 47.56e6; % in Hz
lambda = .5/res; % in m^(-1)
n = (1/(1-(1-Nshots/FOV/lambda)^(1/alpha)));
w = 2*pi*n;
Tea = lambda*w/gamma/Gmax/(alpha+1); % in s
%Tea1 = (lambda*w)/((gamma.*gm)*(alpha+1));
Tes = sqrt(lambda*w^2/SRmax/gamma)/(alpha/2+1); % in s
Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s
 
if Ts2a<Tes
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a);
    Tend = Tea;
else
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
    Tend = Tes;
end
 
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));
dt = Tea*1e-4; % in s
Dt = dt/FOV/abs(k(Tea)-k(Tea-dt)); % in s
t = 0:Dt:Tend; % in s
kt = k(t); % in m^-1
 
DT_GE=4e-6;
tge=0:DT_GE:Tend;
 
    dist = [0,cumsum(abs(kt(2:end)-kt(1:end-1)))];
    kt_x = interp1(t,real(kt),tge,'spline'); 
    kt_y = interp1(t,imag(kt),tge,'spline'); 
    kt_new = complex(kt_x,kt_y);
 
k = kt_new;
end
                
 
 
 
 
 


