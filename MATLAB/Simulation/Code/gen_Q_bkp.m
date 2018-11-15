function [Qpwr_df, Tissue_types,SigmabyRhox,Mass_cell,Mass_corr] = gen_Q(dirname,SAR_type,anatomy)
%% Read input df files and generate Q matrices based on them

disp('Reading deadface files...');
t0_df = cputime;
% Air_seg = df_read(fullfile(dirname,'airSeg.df')); %Tissue density
Tissue_types= df_read(fullfile(dirname,'BodySegmentation.df')); %Tissue density
S = df_read(fullfile(dirname,'material.df')); 
D = size(Tissue_types);
SigmabyRhox = real(S);
Mass_cell =imag(S);
Mass_corr=0;
Mass_air = 1.625e-7 + 1e-9;
clear S S2;

Ex = zeros([size(SigmabyRhox),8]); %dir, vector, tx channels 
Ey=Ex;
Ez=Ex;



%% Lot of hardcoding for now, make it better later

parfor k=1:8 %num_channels
    Ex(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ex.df'])); 
    Ey(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ey.df'])); 
    Ez(:,:,:,k) = df_read(fullfile(dirname,['multix_coil',num2str(k),'_ez.df'])); 

end

clear k  dirname;
fclose('all');
t1_df = cputime - t0_df;
disp(['Deadface files read in ',num2str(t1_df),' seconds']);

switch SAR_type

    case 'global'
           
        %% Based on BodySegmentation.df
              switch anatomy
                  case 'wholebody'
                    R = find(Tissue_types>0); 
                    [X,Y,Z]=ind2sub(size(Tissue_types),R);
        
                  case 'head'
                    R = find(Tissue_types==1); 
                    [X,Y,Z]=ind2sub(size(Tissue_types),R);

                  case 'torso'
                   % TODO

                  case 'extremities'
                   % TODO   
              end

                Qpwr =0;Mass_corr =0;
        
           disp('Performing global Q calculations....');
           dim =2;
           t0_gQ = cputime;
            for r =1:length(X)
             Cr = [X(r),Y(r),Z(r)];
             M = Mass_cell(X(r),Y(r),Z(r));
             
             chk = chkcubair_global(dim,Mass_cell,Mass_air,X(r),Y(r),Z(r));
             
                        if((M > Mass_air) && (chk==1))
                               Qpwr = Qpwr + (M.*(gen_E12ptQ(Ex,Ey,Ez,Cr,SigmabyRhox)));
                               Mass_corr = Mass_corr + M;
                                if(mod(r,1000)==0)
                                    disp (r)
                                end

                        end
            end
            t1_gQ = cputime - t0_gQ;
            disp(['Global Q matrices calculated in ',num2str(t1_gQ),' seconds']);
            disp (r);
            Qpwr_df = Qpwr;
            clear Qpwr;

    case 'local'
            %% Prepare for parfor 
            M= D(1);
            N= D(2);
            P = D(3);
            
        
          ind = zeros(1,M*N*P);    k=0;
          x1=ind;y1=ind;z1=ind;
          disp('Preparing for parfor loop... making indices');
         %% Store indices only which have mass > air, try to make this smarter
         t0_pfp = cputime;
          for x=1:M
              for y=1:N
                  for z=1:P
                             if(Mass_cell(x,y,z) > Mass_air)
                               chk = chkcubair_global(2,Mass_cell,Mass_air,x,y,z);  
                                       if(chk==1)
                                           k=k+1;
                                            ind(k) = sub2ind([M,N,P],x,y,z);
                                            x1(k) =x;y1(k)=y;z1(k)=z;
                                       end
                             end
%                        disp(ind(k));
                  end
              end
          end
          t1_pfp = cputime - t0_pfp;
          disp(['Parfor indices stacked in ',num2str(t1_pfp),' seconds']);
          
          ind = squeeze(ind(ind>0));
          x1 = squeeze(x1(x1>0));
          y1 = squeeze(y1(y1>0));
          z1 = squeeze(z1(z1>0));
          
          Qpwr = zeros(length(ind),8,8); %Remove hardcoding for coils later
          disp('Calculating Qpwr now');
         
          t0_lQ =cputime;
         parfor k=1:length(ind)
                    Qpwr(k,:,:) = gen_E12ptQ(Ex,Ey,Ez,(ind(k)),SigmabyRhox); 
         end
         t1_lQ = cputime - t0_lQ;
         disp(['Local Q matrices calculated in ',num2str(t1_lQ),' seconds']);
    
         clear Ex Ey Ez Tissue_types SigmabyRhox x y z;
         disp('Cleared variables to make space....'); 
            
           
           %% Undo indexing from parfor, overhead - need a smarter implementation 
           Qpwr2 = zeros(M,N,P,8,8);
           tic;
              for k=1:length(ind)
              Qpwr2(x1(k),y1(k),z1(k),:,:) = squeeze(Qpwr(k,:,:));
                    if(mod(k,1e4)==0)
                            disp(k)
                    end
              end
           toc;
           
           clear Qpwr;
           
           %% Average Q over 10g volumes.
           disp('Calculating Mass-averaged local Q matrices...');
            Mdef =0.01; %kg for IEC required mass per an arbitrary volume V
           [Qpwr_df] = squeeze(get_Qavg(Mass_cell,Mdef,Qpwr2,ind,x1,y1,z1));
      
      
                
 end
            
            
           



