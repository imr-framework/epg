%% Calibrate calculations with scanner
S = ones(8)- 2.*eye(8);
% S = hadamard(8);
SAR = zeros(1,length(S));
%%

for k=1:length(S);
    I = conj(squeeze(S(k,:)).');
    SAR(k) = (I'*Q.Qtmf*I)/60;
   
end
 disp(SAR);