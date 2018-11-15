function [SAR,k] = get_SARUI_nocvx(VOPm,x,Bx,By)
% function [SAR,B1server] = get_SARUI(VOPm,x)
ph = 0:pi/4:7*(pi/4);
tx_ph = exp(1i*ph).';
Mc = 0.237125; % for B1 = 5uT

xph = x.*tx_ph;
SAR = zeros(size(VOPm,1),1); 
if(size(VOPm,1)>1)
   for k=1:size(VOPm,1)
        SAR(k) = ((xph)'*squeeze(VOPm(k,:,:))*(xph));
    end
else
    for k=1:size(VOPm,1)
          SAR(k) = ((xph)'*squeeze(VOPm(k,:,:))*(xph));
    end
end
    
figure(101);plot(abs(SAR),'k*-');hold on;drawnow;xlabel('Index');ylabel('Local SAR');
[SAR,k] = max(abs(SAR));
B1server = get_B1server(xph,Bx,By);
SAR = abs((SAR*Mc)/(B1server^2*1e12)); %Normalization by 1e12

