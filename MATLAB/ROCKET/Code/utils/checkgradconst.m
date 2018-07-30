function [channel, dt]= checkgradconst(ktest,dt,SRmax)
gamma =  42576000; %Determined from Pulseq - do not change
ktest = [0 ktest];
G = diff(ktest)/dt/gamma; %T/m
Gx = real(G);
Gy = imag(G);

SRx = diff(Gx)/dt;
SRy = diff(Gy)/dt;

disp([max(abs(SRx)) max(abs(SRy))])

dt1 = diff(Gx)/SRmax;
dt2 = diff(Gy)/SRmax;

dt = max([max(abs(dt1)), max(abs(dt2))]);
disp(dt);
if(max(abs(SRy)) > SRmax)
    channel = 'y';
    [~, ind] = max(abs(SRy));
   
    if(max(abs(SRy)) > SRmax)
        error('Increase number of shots');
    end
    
else
    channel  = 'x';
     [~, ind] = max(abs(SRx));
      if(max(abs(SRx)) > SRmax)
        error('Increase number of shots');
    end
  
end

disp([max(Gx) max(Gy)]);
disp([max(abs(SRx)) max(abs(SRy))]);
