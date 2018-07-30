%% test script

s = linspace(0,100,100);
dt = 1;fact =1;
Gx = (fact./dt).*(diff(squeeze(s)));
Gx = [0 Gx];
stem(s); hold on; stem(Gx);
%%

kx = cumsum(Gx)*dt;

% for l = 1:length(Gx)
%         kx(l) = cumsum(Gx(1:l))*dt(1);     
% end
plot(kx,'k*');