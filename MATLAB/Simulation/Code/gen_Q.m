function [Qpwr_df, Tissue_types,SigmabyRhox,Mass_cell,Mass_corr,Qpwr2] = gen_Q(E,B,M,SAR_type,anatomy)
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


% %% Normalize E fields
% Exm = max(abs(Ex(:)));
% Eym = max(abs(Ey(:)));
% Ezm = max(abs(Ez(:)));
% 
% maxe = max([Exm,Eym,Ezm]);
% Ex = Ex./maxe;
% Ey =Ey./maxe;
% Ez = Ez./maxe;

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
                     R = find(Tissue_types==2); 
                    [X,Y,Z]=ind2sub(size(Tissue_types),R);

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
                                    disp (r);
%                                     figure(101); imagesc(abs(Qpwr)); pause(0.1);
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
            
        
        disp('Preparing for parfor loop... making indices');
         %% Store indices only which have mass > air, try to make this smarter
         ind = Mass_cell > Mass_air;
         [ms]=find(ind);
         disp('Calculating Qpwr now');
         Qpwr = zeros(length(ms),8,8); %Remove hardcoding for coils later

        t0_lQ =cputime;
         
        parfor k=1:length(ms)
                    [m,n,p]=ind2sub(D,ms(k));    
                    Qpwr(k,:,:) = gen_E12ptQ(Ex,Ey,Ez, [m,n,p],SigmabyRhox); 
         end
        
         disp('Creating the 5D Qpwr matrix...');
         Qpwr2 = zeros(M*N*P,8,8);   
         Qpwr2(ms,:,:) = squeeze(Qpwr(1:end,:,:));
         Qpwr2 = reshape(Qpwr2,[M,N,P,8,8]); 

           
         t1_lQ = cputime - t0_lQ;
         disp(['Local Q matrices calculated in ',num2str(t1_lQ),' seconds']);
    
%          clear Ex Ey Ez Tissue_types SigmabyRhox SAR_type Mass_corr anatomy Qpwr t0_df t1_df t0_lQ t1_lQ;
         disp('Cleared variables to make space....'); 

         
         
         
           %% Average Q over 10g volumes.
           disp('Calculating Mass-averaged local Q matrices...');
            Mdef =0.01; %kg for IEC required mass per an arbitrary volume V
           [Qpwr_df] = squeeze(get_Qavg(Mass_cell,Mdef,Qpwr2,ms));
      
      
                
 end
            
            
           



