function radp = get_radkparams(dx,dy,L,dim)

switch dim
    case '2D'
        
            Nx = round(L./dx);
            Ny = round(L/dy);
            dkx  = 1./(Nx.*dx);
            dky = 1./(Ny.*dy);
            kmax = 1./(2.*dx);

            radp.Ns = round(pi.*kmax.*L);
            Nphase = 2*kmax*L;
            radp.Leff = (radp.Ns*2*L/(pi*Nphase));
    case '3D'
           kmax = 1/(2.*dx);
           radp.Ns = round(4*pi*(L.*kmax).^2);
           dtheta = 1/(kmax.*L);
           radp.Ntheta =  ceil(pi/dtheta); % 0 to pi
           radp.Nphi = ceil(radp.Ns/radp.Ntheta); %0 to 2pi
end
        
        
        
        
        
        
        
        
