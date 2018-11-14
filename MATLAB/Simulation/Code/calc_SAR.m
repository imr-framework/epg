function [SAR] = calc_SAR(Q,I)

SAR_temp = zeros(1,size(Q,1));
% n = SAR;
parfor k=1:size(Q,1)
        Qtemp = squeeze(Q(k,:,:));
        SAR_temp(k) = I'*Qtemp*I;
%         n(k) = norm(Qtemp);
end
SAR = max(abs(SAR_temp));