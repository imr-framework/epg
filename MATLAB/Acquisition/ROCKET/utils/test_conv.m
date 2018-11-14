%%
s = zeros(size(x));
s(60:68) = 1;
x = conv(x,s,'same');
X = fftshift(fft(x));
R = fftshift(fft(r2));
H = R.*conj(X);
figure; plot(abs(H));

%%

figure;plot(abs(fft(s)));