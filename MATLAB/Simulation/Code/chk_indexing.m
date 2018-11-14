%% Reading indices based on index file
torso = find(squeeze(Qavg.index(4,:))==1);
x2 =x0_global;
% [theta,rho]=cart2pol(real(x2),imag(x2));
% ph_tx = 0:pi/4:((pi/4)*7);
% theta = theta + ph_tx;
% [x1,y1]= pol2cart(theta,rho);
% x2 = complex(x1,y1);

ph_tx = 0:pi/4:((pi/4)*7);
ph = exp(1i*ph_tx);
x2 = x2.*ph.';

SAR=zeros(1,length(torso));
for k=1:length(torso)
    ind = torso(k);
    SARt =x2'*squeeze(Qavg.avg(ind,:,:))*x2;
    SAR(k) = sum(abs(SARt(:)));
end

SAR = abs(max(SAR));
disp(SAR);




