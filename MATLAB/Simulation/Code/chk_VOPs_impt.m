%%
A = Exp.Asys;
Ainv = pinv(abs(A));


%%
for k=1:72
    VOP1 = squeeze(VOPm(k,:,:));
    M = (Ainv)'*VOP1*(Ainv);
    imagesc(abs(M));caxis([0 0.25])
    pause;
end

