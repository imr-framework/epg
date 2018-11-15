% function [SAR] = get_SAR(VOPm,x)
function [SAR,B1server] = get_SARUI(VOPm,x)
ph = 0:pi/4:7*(pi/4);
tx_ph = exp(1i*ph).';
Mc = 0.237125; % for B1 = 5uT


% SAR = zeros(size(VOPm,1),1); to maintain consistency b/w MATLAB & cvx
% Can make SAR an expression but is problematic.
if(size(VOPm,1)>1)
    
    for k=1:size(VOPm,1)
    %     SAR(k) = (x'*squeeze(VOPm(k,:,:))*x);
    SAR(k) = ((x.*tx_ph)'*squeeze(VOPm(k,:,:))*(x.*tx_ph));
    end
else
    for k=1:size(VOPm,1)
    %     SAR(k) = (x'*squeeze(VOPm(k,:,:))*x);
    SAR(k) = ((x.*tx_ph)'*squeeze(VOPm(k,:,:))*(x.*tx_ph));
    end
end


% %% Change this back to original after this experiment
% SAR = ((x.*tx_ph)'*squeeze(VOPm(1,:,:))*(x.*tx_ph));
[SAR] = max(SAR);

B1server = get_B1server(x); 
SAR = (SAR*Mc)/(B1server^2*1e12); %Normalization by 1e12

