function Q = chkcubair(dim,Mass_cell,Mass_air,x,y,z,Q)
%% Check if cube's face is in air
               
                xface0 = squeeze(Mass_cell(x-dim,y-dim:y+dim,z-dim:z+dim));
                yface0 = squeeze(Mass_cell(x-dim:x+dim,y-dim,z-dim:z+dim));
                zface0 = squeeze(Mass_cell(x-dim:x+dim,y-dim:y+dim,z-dim));
                
                
                xface1 = squeeze(Mass_cell(x+dim,y-dim:y+dim,z-dim:z+dim));
                yface1 = squeeze(Mass_cell(x-dim:x+dim,y+dim,z-dim:z+dim));
                zface1 = squeeze(Mass_cell(x-dim:x+dim,y-dim:y+dim,z+dim));
                
                
                
                
                
               if((max(xface0(:)) <= Mass_air)||(max(yface0(:)) <= Mass_air)||(max(zface0(:)) <= Mass_air) ||...
                      (max(xface1(:)) <= Mass_air)||(max(yface1(:)) <= Mass_air)||(max(zface1(:)) <= Mass_air))  
                    disp('Invalid cell to calculate SAR');
               end
                
%% All those cells which used this cell have to be set to max value
% All cells in the range of 2*dim + 1 of this cell would have used this

Qsump =0;

for i =    x-2*dim:x+2*dim
    for j =  y-2*dim:y+2*dim
        for k = z-2*dim:z+2*dim
   
        S = squeeze(Q(i,j,k,:));
        Qsumc = sum(abs(S(:)));
        
        if(Qsumc > Qsump)
            Qsump = Qsumc;
            maxi=i;
            maxj=j;
            maxk=k;
        end
            
            
            
            
        end
    end
end

%% Set those cells which used this cell to average to the peak value
% Q(x-2*dim:x+2*dim,y-2*dim:y+2*dim,z-2*dim:z+2*dim) = Q(maxi,maxj,maxk,:,:);          
%TODO - vectorize this for loop 
for i = x-2*dim:x+2*dim
    for j = y-2*dim:y+2*dim
        for k = z-2*dim:z+2*dim
   
        if(Mass_cell(i,j,k) > Mass_air)
            Q(i,j,k,:,:) = Q(maxi,maxj,maxk,:,:);
        end
            
            
            
        end
    end
end
               
                
%% if one face is is air, then set the value of the cells that used this cell to the maximum 
Q(x,y,z,:,:) = -1;         