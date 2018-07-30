function kx = checkfix_k(kx,delk_max)

ind = diff(kx)>delk_max;
if(sum(ind) >0)
    for k = 1:length(ind)
        
    end
end