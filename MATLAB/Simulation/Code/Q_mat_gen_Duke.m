%% Read the required data in
% phantom data has e fields files named as_Ex while the VHM has _ex

% Author: Sairam Geethanath, Ph.D.
% 01/2012
% phantom data has e fields files named as_Ex while the VHM has _ex
addpath(genpath('D:\sar_code_VOP_SG'));
close all;
SAR_type = 'Global';
% if(parpool('local') == 0)
%   parpool local 8;
% end
Mass_body = 71.5;%Duke


%%
t0 = cputime;
switch SAR_type 

    case 'Global'
            clc;disp('SAR type: GLOBAL');
           
           
            %% Read from df files and check local computation, using index - whole body
            dirname = uigetdir('');
            clc;disp('Q - Whole body calculation started ....');
            [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body]  = gen_Q_custom(dirname,'global','wholebody');%SbR stands for SigmabyRho
            Qavg_tm = Qavg_df./Mass_body;
            figure(1);subplot(122);imagesc(abs((Qavg_tm)));colorbar;title('Implemented - Mass normalized BODY');
            figure(3);imagesc(abs(abs(Q.Qtmf) - abs(Qavg_tm)));colorbar;title('Difference BODY');
            figure(4);imagesc(abs(abs(Q.Qtmf)./abs(Qavg_tm)));colorbar;title('Ratio BODY');

%             %% Head
%             disp('Q -Head calculation started ....');
%             [Qavg_df,~,~,~,Mass_head] = gen_Q(dirname,'global','head');%SbR stands for SigmabyRho
%             Qavg_hm = Qavg_df./Mass_head;
%             figure(2);subplot(122);imagesc(abs(Qavg_hm));colorbar;title('Implemented - Mass normalized HEAD');
%             figure(5);imagesc(abs(abs(Q.Qhmf) - abs(Qavg_hm)));colorbar;title('Difference HEAD');
%             figure(6);imagesc(abs(abs(Q.Qhmf)./  abs(Qavg_hm)));colorbar;title('Ratio HEAD');
d
            
           
            
    
    case 'Local'
            disp('SAR type: LOCAL');
    
            %%
            dirname = uigetdir('');   
            disp('Starting calculation of local Q matrices.....');
            t0_local = cputime;
                 [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body] = gen_Q_custom(dirname,'local','wholebody');%SbR stands for SigmabyRho
            t1_local = cputime - t0_local;
            clc;disp(['Done calculating local Q-matrices in ',num2str(t1_local),' seconds']);
           
            
                %% Save
                Qavg.imp= Qavg_df;
                save('LocalQ','Qavg','-v7.3'); %Could do uisave as well, but for building VOPs....


end
matlabpool close;
t1 = cputime - t0;
disp(['Q matrices calculated in ',num2str(t1),' seconds']);



