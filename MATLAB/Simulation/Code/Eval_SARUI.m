%% Evaluation of the implementations of B1server and SAR UI calculations
% Working only on philips based Q matrices
%% Load required files
clear all;
first_run =1;
if(first_run)
    Seg_areas = df_read('BodySegmentation_human'); %1 - brain, 2 - torso, 3 - extremities
    [Filename,Pathname]=uigetfile('*.mat','Load Global Q matrix'); %choose head for relevant simulation
    load(fullfile(Pathname,Filename));

    [Filename,Pathname]=uigetfile('*.mat','Load local Q matrix'); %choose head for relevant simulation
    load(fullfile(Pathname,Filename));
    
    Sarlog_row_header = {'Whole body','Exposed Mass','Global Head','Local Head','Local Torso','Local Extremeties','Bplus'};
    xlswrite('Sarlog_predicted.xls',Sarlog_row_header);
    row =1;
end



%% Pick shims and perform exhaustive search
% x = 1                                       .5.*ones(8,1);
x = [0;0;1;1;1;1;0;0];

%%  Global search
%Global total body
Q(1,:,:) = Qavg_philips.Qtmf;
[SAR_tbody,B1server] = get_SARUI(Q,x);

%Global head
clear Q;
Q(1,:,:) = Qavg_philips.Qhmf;
[SAR_thead, B1server] = get_SARUI(Q,x);


%Global exposed body
clear Q;
Q(1,:,:) = Qavg_philips.Qemf;
[SAR_texposed,B1server] = get_SARUI(Q,x);



%% Local search
clear Q;
index = squeeze(Qavg.index(4,:));

% Local head 
ind_head = find(index==1);
Q = squeeze(Qavg.avg(ind_head,:,:));
[SAR_lhead, B1server] = get_SARUI(Q,x);


% Local torso 
ind_torso = find(index==2);
Q = squeeze(Qavg.avg(ind_torso,:,:));
[SAR_ltorso,B1server] = get_SARUI(Q,x);


% Local torso 
ind_extremeties = find(index==3);
Q = squeeze(Qavg.avg(ind_extremeties,:,:));
[SAR_lextremeties,B1server] = get_SARUI(Q,x);


Sarlog_row = {abs(SAR_tbody),abs(SAR_texposed),abs(SAR_thead),abs(SAR_lhead),abs(SAR_ltorso),abs(SAR_lextremeties),abs(B1server)}; %Missing shim info
row = row + 1;
xlswrite('Sarlog_predicted.xls',Sarlog_row,1,['A',num2str(row)]);

%%
writeshims_r253_sjm(x,'nonorm','name','b1test_ rfshims');