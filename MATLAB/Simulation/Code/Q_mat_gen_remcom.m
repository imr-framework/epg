%% Read the required data in
% Author: Sairam Geethanath, Ph.D.
% 01/2018
% phantom data has e fields files named as_Ex while the VHM has _ex
addpath(genpath('D:\sar_code_VOP_SG'));
close all;
SAR_type = 'Global';

%% Load the relevant REMCOM data file from Jinfeng

dirname = uigetdir('','Pick the directory that has E-fields, B-fields and M');
[E,B,M] = prep_data(dirname);
%%
t0 = cputime;
switch SAR_type 

    case 'Global'
            clc;disp('SAR type: GLOBAL');
           
            %% Read from df files and check local computation, using index - whole body
            clc;disp('Q - Whole body calculation started ....');
            [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body]  = gen_Q_remcom(E,M,'global','wholebody');%SbR stands for SigmabyRho
            Qavg_tm = Qavg_df./Mass_body;
            figure(1);imagesc(abs((Qavg_tm)));colorbar;title('Implemented - Mass normalized BODY');
                
%             figure(1);subplot(122);imagesc(abs((Qavg_tm)));colorbar;title('Implemented - Mass normalized BODY');
%             figure(3);imagesc(abs(abs(Q.Qtmf) - abs(Qavg_tm)));colorbar;title('Difference BODY');
%             figure(4);imagesc(abs(abs(Q.Qtmf)./abs(Qavg_tm)));colorbar;title('Ratio BODY');

            %% Head - Uncomment once the body work starts
%             disp('Q -Head calculation started ....');
%             [Qavg_df,~,~,~,Mass_head] = gen_Q(dirname,'global','head');%SbR stands for SigmabyRho
%             Qavg_hm = Qavg_df./Mass_head;
%             figure(2);subplot(122);imagesc(abs(Qavg_hm));colorbar;title('Implemented - Mass normalized HEAD');
%             figure(5);imagesc(abs(abs(Q.Qhmf) - abs(Qavg_hm)));colorbar;title('Difference HEAD');
%             figure(6);imagesc(abs(abs(Q.Qhmf)./abs(Qavg_hm)));colorbar;title('Ratio HEAD');
    
            
            %% Make structure for writing it to a file
            Q.Qtmf = Qavg_tm;
            save('GlobalQ','Q','-v7.3');
%             Q.Qhmf = Qavg_hm;
%             Q.Qemf = Q.Qemf; %TODO have not implemented this at all, need to figure out
%             [status] = write_qmat(Q,'Global');
    
    case 'Local'
              [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body]  = gen_Q_remcom(E,M,'local','wholebody');%SbR stands for SigmabyRho
              Qavg_tm = Qavg_df./Mass_body;
%                figure(1);imagesc(abs((Qavg_tm)));colorbar;title('Implemented - Mass normalized BODY');

%             disp('Starting calculation of local Q matrices.....');
%             t0_local = cputime;
%             [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body] = gen_Q(dirname,'local','wholebody');%SbR stands for SigmabyRho
%             t1_local = cputime - t0_local;
%             clc;disp(['Done calculating local Q-matrices in ',num2str(t1_local),' seconds']);
%            
%             Qpt_design = squeeze(Qavg_df(61,50,42,:,:));
%             figure;imagesc(abs((Qpt_design)));
%             figure;imagesc(abs(Qpt)./abs((Qpt_design)));
%             figure;imagesc(abs(abs(Qpt) -abs((Qpt_design))));
%             %%
%             [RMSE_map, Qtri_map, Qimp_map,NRMSE_map] = get_RMSE(Qavg,Qavg_df,Mass_cell);
%             
%             %% Global Display of local Q
%                 figure;imshow(squeeze(Qtri_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 2]);
%                 figure;imshow(squeeze(Qimp_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 2]);
%                 diff_map = abs(Qtri_map - Qimp_map);
%                 figure;imshow(squeeze(diff_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 0.1])

                %% Save
                Qavg.imp= Qavg_df;
              
                 save('LocalQ','Qavg','-v7.3');


end
% parpool close;
t1 = cputime - t0;
disp(['Q matrices calculated in ',num2str(t1),' seconds']);



