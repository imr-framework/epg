function B1 = get_gapped_DANTE(A0,beta,myu, Nsamp,L,system)

Tp = 2*Nsamp.*system.gradRasterTime; %50% of the time is pulse and the other half is the adc = 2560 us
t = -Tp/2:system.gradRasterTime:Tp/2;
[x,A,w1] = make_HSI(A0,beta,myu,t,1); %AFP
r = zeros(size(t));
r(128:129) =1;

xp = conv(x,r,'same');
xgap = zeros(size(xp));
xgap(1:2*L:end)  = xp(1:2*L:end);
xgapp = conv(xgap, r,'same');



