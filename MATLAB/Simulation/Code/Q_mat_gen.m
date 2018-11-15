%% Read the required data in
% Author: Sairam Geethanath, Ph.D.
% 01/2012
% phantom data has e fields files named as_Ex while the VHM has _ex
addpath(genpath('D:\sar_code_VOP_SG'));
close all;
SAR_type = 'Global';
% if(parpool('local') == 0)
%   parpool local 8;
% end

%%
t0 = cputime;
switch SAR_type 

    case 'Global'
            clc;disp('SAR type: GLOBAL');
            %% Get Philips Q values
            [Q] = read_qmat('globalQMatrix_human');%             Qtmf = Q.Qtmf;         Qhmf = Q.Qhmf;             Qemf = Q.Qemf;
            figure(1);subplot(121);imagesc(abs(Q.Qtmf));colorbar;title('Whole body');
            figure(2);subplot(121);imagesc(abs(Q.Qhmf));colorbar;title('Head');
           
            %% Read from df files and check local computation, using index - whole body
            dirname = uigetdir('');
            clc;disp('Q - Whole body calculation started ....');
            [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body]  = gen_Q_bkp(dirname,'global','wholebody');%SbR stands for SigmabyRho
            Qavg_tm = Qavg_df./Mass_body;
            figure(1);subplot(122);imagesc(abs((Qavg_tm)));colorbar;title('Implemented - Mass normalized BODY');
            figure(3);imagesc(abs(abs(Q.Qtmf) - abs(Qavg_tm)));colorbar;title('Difference BODY');
            figure(4);imagesc(abs(abs(Q.Qtmf)./abs(Qavg_tm)));colorbar;title('Ratio BODY');

            %% Head
            disp('Q -Head calculation started ....');
            [Qavg_df,~,~,~,Mass_head] = gen_Q(dirname,'global','head');%SbR stands for SigmabyRho
            Qavg_hm = Qavg_df./Mass_head;
            figure(2);subplot(122);imagesc(abs(Qavg_hm));colorbar;title('Implemented - Mass normalized HEAD');
            figure(5);imagesc(abs(abs(Q.Qhmf) - abs(Qavg_hm)));colorbar;title('Difference HEAD');
            figure(6);imagesc(abs(abs(Q.Qhmf)./  abs(Qavg_hm)));colorbar;title('Ratio HEAD');

            
             %% Head
%             disp('Q - torso calculation started ....');
%             [Qavg_df,~,~,~,Mass_torso] = gen_Q(dirname,'global','torso');%SbR stands for SigmabyRho
%             Qavg_hm = Qavg_df./Mass_torso;
%             figure(2);subplot(122);imagesc(abs(Qavg_hm));colorbar;title('Implemented - Mass normalized HEAD');
%             figure(5);imagesc(abs(abs(Q.Qhmf) - abs(Qavg_hm)));colorbar;title('Difference HEAD');
%             figure(6);imagesc(abs(abs(Q.Qhmf)./  abs(Qavg_hm)));colorbar;title('Ratio HEAD');

            
            
            %% Make structure for writing it to a file
            Q.Qtmf = Qavg_tm;
            Q.Qhmf = Qavg_hm;
            Q.Qemf = Q.Qemf; %TODO have not implemented this at all, need to figure out
            [status] = write_qmat(Q,'Global');
    
    case 'Local'
            disp('SAR type: LOCAL');
            [Qavg] = read_qmat('avgQMatrix');
            %%
            S = squeeze(Qavg.index(1:3,:));
            
            x =find(S(1,:)== 60);
            y =find(S(2,:)==49);
            z =find(S(3,:)==42);
            
            s = intersect(x,y);
            s = intersect(s,z);
            
            Qpt = squeeze(Qavg.avg(s,:,:));
            figure;imagesc(abs(Qpt));

            dirname = uigetdir('');   
            disp('Starting calculation of local Q matrices.....');
            t0_local = cputime;
                 [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body] = gen_Q(dirname,'local','wholebody');%SbR stands for SigmabyRho
            t1_local = cputime - t0_local;
            clc;disp(['Done calculating local Q-matrices in ',num2str(t1_local),' seconds']);
           
            Qpt_design = squeeze(Qavg_df(61,50,42,:,:));
            figure;imagesc(abs((Qpt_design)));
            figure;imagesc(abs(Qpt)./abs((Qpt_design)));
            figure;imagesc(abs(abs(Qpt) -abs((Qpt_design))));
            %%
            [RMSE_map, Qtri_map, Qimp_map,NRMSE_map] = get_RMSE(Qavg,Qavg_df,Mass_cell);
            
            %% Global Display of local Q
                figure;imshow(squeeze(Qtri_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 2]);
                figure;imshow(squeeze(Qimp_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 2]);
                diff_map = abs(Qtri_map - Qimp_map);
                figure;imshow(squeeze(diff_map(:,:,193)),'InitialMagnification',200);colorbar;colormap(jet);caxis([0 0.1])

                %% Save
                Qavg.imp= Qavg_df;
                 save('LocalQ','Qavg','-v7.3');


end
parpool close;
t1 = cputime - t0;
disp(['Q matrices calculated in ',num2str(t1),' seconds']);



