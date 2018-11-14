function [E,B,M] = prep_data(dirname)

%% This file prepares data for the SAR calculations code at ICL
disp('Reading MATLAB files  ...');
t0_df = cputime;
pname=dirname;
%% Read E fields
% [fname,pname]=uigetfile('*.mat','Pick the E fields in Mat format');
% fname = 'Efields.mat';
% load(fullfile(pname,fname));
% E = permute(CEFT, [2 3 4 5 1]);
[E,B] = load_EMfields(pname);
%% Read B fields
% [fname,pname]=uigetfile('*.mat','Pick the B fields in Mat format');
% B = load_Bfields(pname);
% fname = 'Bfields.mat';
% load(fullfile(pname,fname));
% B = permute(B1PT, [2 3 4 5 1]);

%% Read Material properties - this needs to be better defined in the future - TODO Jinfeng
% [fname,pname]=uigetfile('*.mat','Pick the Mask in 4D');
fname ='Materials.mat';
load(fullfile(pname,fname));
M.Tissue_types = Mask_Phantom(:,:,:,1).*Mask_Phantom(:,:,:,2).*Mask_Phantom(:,:,:,3); %This needs to be better coded for later
Volume_cell = (2*1e-3).^3; %This is 2mm x 2mm x 2mm but needs to be from a file later
Rhox = 1000; %kg/m^3; %This needs to be from a file
Sigma = 0.64;%S/m - This needs to be from a file
M.SigmabyRhox = Sigma/Rhox;  %0.64 S/m/kg/m^3 over This needs to be from a file later
M.Mass_air = 1.625e-7 + 1e-9;
M.Mass_cell = Volume_cell.*Rhox;%mass = volume.*density

