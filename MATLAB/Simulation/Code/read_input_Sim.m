function [Qpwr, Tissue_types,SigmabyRhox,Mass_cell,Mass_corr] = read_input_Sim(dirname,SAR_type,index,anatomy)
%% Read input df files and generate Q matrices based on them

Tissue_types= df_read(fullfile(dirname,'BodySegmentation.df')); %Tissue density
% Tissue_types= df_read(fullfile('BodySegmentation_human.df')); %Tissue density
S = df_read(fullfile(dirname,'material.df')); 
S2 = df_read(fullfile(dirname,'material_yz.df')); 

SigmabyRhox = real(S);
Mass_cell =imag(S);
SigmabyRhoy = real(S2);
SigmabyRhoz = imag(S2);
Mass_corr=0;

% %% Employ index file to delineate head voxels.
% R = find(squeeze(index(4,:))>0); %1 is head
% X=squeeze(index(1,R));
% Y=squeeze(index(2,R));
% Z=squeeze(index(3,R));

clear S S2;

Ex = zeros([size(SigmabyRhox),8]); %dir, vector, tx channels 
Ey=Ex;
Ez=Ex;



%% Lot of hardcoding for now, make it better later

for k=1:8 %num_channels
    Ex(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ex.df'])); 
    Ey(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ey.df'])); 
    Ez(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ez.df'])); 
end
clear k  dirname;
fclose('all');

% Qpwr = zeros(8,8,length(X),length(Y),length(Z));

switch SAR_type

    case 'global'
           
        %% Based on BodySegmentation.df
              switch anatomy
                  case 'wholebody'
%                     R = find(squeeze(index(4,:))>0); 
%                     X=squeeze(index(1,R));
%                     Y=squeeze(index(2,R));
%                     Z=squeeze(index(3,R));

                    R = find(Tissue_types>0); 
                    [X,Y,Z]=ind2sub(size(Tissue_types),R);
        
                  case 'head'
%                     R = find(squeeze(index(4,:))==1); 
%                     X=squeeze(index(1,R));
%                     Y=squeeze(index(2,R));
%                     Z=squeeze(index(3,R));

                    R = find(Tissue_types==1); 
                    [X,Y,Z]=ind2sub(size(Tissue_types),R);

                  case 'torso'
%                     R = find(squeeze(index(4,:))==2); 
%                     X=squeeze(index(1,R));
%                     Y=squeeze(index(2,R));
%                     Z=squeeze(index(3,R));
                  case 'extremities'
%                     R = find(squeeze(index(4,:))==3);
%                     X=squeeze(index(1,R));
%                     Y=squeeze(index(2,R));
%                     Z=squeeze(index(3,R));
              end

        
        
        
        Qpwr =0;Mass_corr =0;
        
        
        % %% Employ index file to delineate head voxels.
       
        
        
        
           disp('Performing global Q calculations....');
            for r =1:length(X)
             Cr = [X(r),Y(r),Z(r)];
             M = Mass_cell(X(r),Y(r),Z(r));
             Qpwr = Qpwr + (M.*(gen_E12ptQ(Ex,Ey,Ez,Cr,SigmabyRhox,SigmabyRhoy,SigmabyRhoz)));
             Mass_corr = Mass_corr + M;
                if(mod(r,1000)==0)
                 disp (r)
                end
            end
            disp (r);
    case 'local'
             Qpwr = zeros(60,60,60,8,8); %need to compute atleast 125 values for 1 location
%             for r =1:length(X)
%              Cr = [X(r),Y(r),Z(r)];
              
              pt = [30,30,18];
              
              for x=1:5
                  for y=1:5
                      for z=1:5
                        Cr = [pt(1)-3+x, pt(2)-3+y,pt(3)-3+z];
                        Qpwr(pt(1)-3+x,pt(2)-3+y,pt(3)-3+z,:,:) = (gen_E12ptQ(Ex,Ey,Ez,Cr,SigmabyRhox,SigmabyRhox,SigmabyRhox));

                      end
                  end
              end
              

          %% Average Q over 10g volumes.
           Mdef =0.01; %kg for IEC required mass per an arbitrary volume V
           clear Qpwr;
           [Qpwr] = squeeze(get_Qavg(Mass_cell,Mdef,Qpwr));
      
        

           
           
           
           
           
           
           
           
           
        
%         %%
           
            
                
 end
            
            
           


% function M = gen_M12pt(Mass_cell,R)
