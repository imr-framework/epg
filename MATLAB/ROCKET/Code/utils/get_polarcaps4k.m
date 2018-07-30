function k = get_polarcaps4k(k,gamma,system,kind)

kx_target = k(1,kind); 
ky_target = k(2,kind);
kz_target = k(3,kind);

%% Find an archimedian spiral that conforms to system and generatees the required kmax
klim_target = sqrt(kx_target.^2 + ky_target.^2);
alpha =1;
Nx = 
[ktraj,G,lambda]= vds2D_pulseq_v1(fov,Nx,Nshots,alpha,system);

