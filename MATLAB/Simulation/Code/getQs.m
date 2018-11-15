function [Qt]= getQs(Q,Mass_cell,dim,x,y,z,f,Mdef)

Qt=0;
 
for i = -dim:dim
    for j=-dim:dim
        for k = -dim:dim
   
                x0 = x - dim;
                y0 = y - dim;
                z0 = z -dim;
                
                M = getm(Mass_cell,x0,y0,z0);
                            status = 0;
							
                            if (abs(i) == dim)
							   status = status +1;
                            end
                            
							if (abs(j) == dim)
							   status = status +1;
                            end
                            
							if (abs(k) == dim)
							   status = status +1;
                            end
                            
                           
                            switch status
                                
                                case 0 %Core
                                    Qt = Qt + (squeeze(Q(x0,y0,z0,:,:)).*M);
                             
                                case 1 %Face
                                    Qt = Qt + (squeeze(Q(x0,y0,z0,:,:)).*M.*f);
                                   
                                case 2 %Edge 
                                    Qt = Qt + (squeeze(Q(x0,y0,z0,:,:)).*M.*f^2);
                                    
                                case 3 %Corner 
                                    Qt = Qt + (squeeze(Q(x0,y0,z0,:,:)).*M.*f^3);
                              
                            end
        
        end
    end
end

Qt = Qt./(Mdef);


function M = getm(Mass_cell,x,y,z)

 M = (1/12).*((3*Mass_cell(x,y,z)) + (2*Mass_cell(x+1,y,z)) + (2*Mass_cell(x,y+1,z)) + (2*Mass_cell(x,y,z+1)) + ...
                    Mass_cell(x+1,y+1,z) + Mass_cell(x+1,y,z+1) + Mass_cell(x,y+1,z+1));