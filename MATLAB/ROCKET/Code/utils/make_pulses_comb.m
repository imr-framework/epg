function [tplist, rfrdlist] = make_pulses_comb()

% Tfull = 30; %s
% Npulses = 128; %Fixed for Siemens 

tpmin = 10e-6;%s
rfrdmin = 10e-6;%s
rfraster = 1e-6; %s default scanner spec

%%
tplist(1) = tpmin;
rfrdlist(1) = rfrdmin;
ct = 1;
for p=1:8
   tplist(ct) = tpmin - rfraster;
  for q = 1: 16
              ct = ct+1;
              tplist(ct) = tplist(ct-1) + rfraster;
              rfrdlist(ct) = rfrdlist(ct-1);
%               disp([tplist(ct)  rfrdlist(ct)])
    
      
  end 
    rfrdlist(ct) = rfrdlist(ct -1) + rfraster;
end

tplist = tplist(1:128);
rfrdlist = rfrdlist(1:128);
