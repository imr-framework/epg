function [B1,A,w1] = make_HSI(A0,beta,myu,t,disp)
A = A0.*sech(beta.*t);
w1 = -myu.*beta.*tanh(beta.*t);
B1 = A.*exp(-1i*w1.*t);
% a = atan(A./w1);
% adiabatic_cond = diff(a)./diff(t);


%%
if(disp)
figure(211); subplot(511); stem(t,A,'k'); xlabel('time(ms)'); ylabel('AM envelope'); hold all;
subplot(512); plot(t,w1,'k');xlabel('time(ms)'); ylabel('FM envelope');hold all;
subplot(513);plot(t,real(B1)); hold on; plot(t, imag(B1)); xlabel('time(ms)'); ylabel ('B1'); legend('B1x','B1y');
subplot(514);plot(real(B1), imag(B1)); xlabel('B1x'); ylabel('B1y');
subplot(515); plot(A,w1); xlabel(''); ylabel('myuB');hold on;
% figure(2); plot(a);
end