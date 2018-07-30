function  k = fix_ktraj(k,indk,system,gamma)

dbstop if error;
fact=2;
count =0;
if(indk <9)
    indk=9;
end
while(indk>0)
    if(indk <9)
        indk=9;
    end
    kn=[];count = count +1;
    kx_append = interp(squeeze(k(1,1:indk+5)),fact); kn(1,:) = [kx_append k(1,indk+6:end)];
    ky_append = interp(squeeze(k(2,1:indk+5)),fact); kn(2,:) = [ky_append k(2,indk+6:end)];
    kz_append = interp(squeeze(k(3,1:indk+5)),fact); kn(3,:) = [kz_append k(3,indk+6:end)];
    indk = get_ktraj2SR(kn, gamma, system, 0);
    k = kn;
    if(count > 20)
        error('This is getting computationally intractable');
    end
end