function VOP_map = map_vops(cluster,ind_sorta,ind)

cluster = squeeze(cluster(cluster>0));
VOP_map = zeros(60^3,1);
label=1;
st =0;
for k=1:length(cluster)
    num_vox = cluster(k);
    
    
    VOP(ind(ind_sorta(st+1:num_vox+st)) = label;
    
    
end