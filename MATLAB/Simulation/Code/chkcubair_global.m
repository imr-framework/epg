function chk = chkcubair_global(dim,Mass_cell,Mass_air,x,y,z)
%% Check if cube's face is in air
               
                xface0 = squeeze(Mass_cell(x-dim,y-dim:y+dim,z-dim:z+dim));
                yface0 = squeeze(Mass_cell(x-dim:x+dim,y-dim,z-dim:z+dim));
                zface0 = squeeze(Mass_cell(x-dim:x+dim,y-dim:y+dim,z-dim));
                
                
                xface1 = squeeze(Mass_cell(x+dim,y-dim:y+dim,z-dim:z+dim));
                yface1 = squeeze(Mass_cell(x-dim:x+dim,y+dim,z-dim:z+dim));
                zface1 = squeeze(Mass_cell(x-dim:x+dim,y-dim:y+dim,z+dim));
                
                chk=1;
                
                
                
               if((max(xface0(:)) <= Mass_air)||(max(yface0(:)) <= Mass_air)||(max(zface0(:)) <= Mass_air) ||...
                      (max(xface1(:)) <= Mass_air)||(max(yface1(:)) <= Mass_air)||(max(zface1(:)) <= Mass_air))  
                       chk=0;
%                       disp('Invalid cell to calculate SAR');
               end
                
