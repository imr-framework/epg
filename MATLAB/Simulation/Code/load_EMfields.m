function [E,B1P] = load_EMfields(pname)

files = dir(pname);
nch = 8; %Number of transmit channels
temp = load(fullfile(pname,files(4).name));

%% Initiate variables for storing
E = zeros(size(temp.TotalField_E_X,1), size(temp.TotalField_E_X,2), size(temp.TotalField_E_X,3),3, nch);
B1P = E;

for ch = 1: nch % 1 B1P file, 1 E file for the same channel
    % B1p
    B1pname = ['Ch0',num2str(ch),'_B1P.mat'];
    load(fullfile(pname,B1pname));
    
    B1P(:,:,:,1,ch) = TotalField_RotatingBPlus_X;
    B1P(:,:,:,2,ch) = TotalField_RotatingBPlus_Y;
    B1P(:,:,:,3,ch) = TotalField_RotatingBPlus_Z;
    
    % E fields 
    Ename = ['Ch0',num2str(ch),'_E.mat'];
    load(fullfile(pname,Ename));
    E(:,:,:,1,ch) = TotalField_E_X;
    E(:,:,:,2,ch) = TotalField_E_Y;
    E(:,:,:,3,ch) = TotalField_E_Z;
 
end

%% Permute dimensions to undo Remcom encoding of Z,Y,X to X,Y,Z
B1P = permute(B1P, [3 2 1 4 5]);
E = permute(E,[3 2 1 4 5]);


