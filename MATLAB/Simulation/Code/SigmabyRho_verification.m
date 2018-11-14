%% Read the required data in
% phantom data has e fields files named as_Ex while the VHM has _ex

close all;
SAR_type = 'Local';
if(matlabpool('size') == 0)
  matlabpool local 4;
end
%%
dirname = uigetdir('','Pick normed fields directory');
      t0_local = cputime;
               [Qavg_df,Tissue_types,SbRx,Mass_cell,Mass_body,Qpwr2] = gen_Q(dirname,'local','wholebody');%SbR stands for SigmabyRho
       t1_local = cputime - t0_local;

[Filename,Pathname]=uigetfile('*.mat','Pick the local Q matrices');
load(fullfile(Pathname,Filename));

%%



