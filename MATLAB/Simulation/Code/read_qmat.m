function [Q] = read_qmat(fname)

%% Read global matrix
% fid = fopen('globalQMatrix.tri','r','l');
fname1 = [fname,'.tri'];
fid = fopen(fname1,'r','l');
header = fread(fid,'uint32');
header = squeeze(header(1:3));
fclose(fid);

fid = fopen(fname1,'r','l');
data = fread(fid,'float');
data = squeeze(data(4:end));
fclose(fid);


%% Header goes - NrofCoils, NrofCoil, NrofCells (4 x3 =12)
% NrofCoils = data1(1);
NrofCoils = header(2);
NrofCells = header(3);
disp('Reading Q matrices....');
switch NrofCells
    
    case 3 %Global SAR q-matrices.
   %% Exploit hardcoding here because these matrix sizes will not change
        R = data(1:2:end);
        I = data(2:2:end);
        
        Qavg = complex(R,I);
        
        %presence of only three cells.
        Qtm = Qavg(1:(end/3));%Qtotalmass (i.e. corresponding to 104.874706 kgs)
        Qhm = Qavg((end/3)+1:(2*end/3));%Qheadmass (i.e. corresponding to 6.503450 kgs)
        Qem = Qavg((2*end/3)+1:(end));%Qexposedmass (i.e. corresponding to 36.666021 kgs)
        
        
        %% Put them as upper triangular matrices.
        Qtmfu = triu(ones(NrofCoils),0); 
        Qtmfu(Qtmfu==1) = Qtm;
        Qtmf  = (triu(Qtmfu,1))' + Qtmfu; %adding upper and lower matrices to form the full matrix, Hermetian.
                                
        Qhmfu = triu(ones(NrofCoils),0); 
        Qhmfu(Qhmfu==1) = Qhm;
        Qhmf  = (triu(Qhmfu,1))' + Qhmfu;
              
        Qemfu = triu(ones(NrofCoils),0); 
        Qemfu(Qemfu==1) = Qem;
        Qemf  = (triu(Qemfu,1))' + Qemfu;
        
        Q.Qtmf = Qtmf;
        Q.Qhmf = Qhmf;
        Q.Qemf = Qemf;

     
        
        
        
    otherwise
     

        
  
        
        
        %% Read associated index file - better to have 2 seperate matrices
        % than one big one.
        disp('Reading index file now....');
        fname = [fname,'.index'];
            fid = fopen(fname,'r','l');
            data_ind = fread(fid,'uint32');
            data_ind = squeeze(data_ind(1:4));
         fclose(fid);

            dimx = data_ind(1);
            dimy =data_ind(2);
            dimz = data_ind(3);
            NrOfSarCells = data_ind(4);
%            NrOfSarCells = (length(data_ind)-8)/4;
             
            fid = fopen(fname,'r','l');
            data_ind = fread(fid,'uint16');
            fclose(fid);


           index = zeros(4,NrOfSarCells);
           index(1,:) = data_ind(9:4:end); %x       
           index(2,:) = data_ind(10:4:end); %y
           index(3,:) = data_ind(11:4:end); %z
           index(4,:) = data_ind(12:4:end); %label
        
        
         

        R = data(1:2:end);
        I = data(2:2:end);
       NrofCells= NrOfSarCells;
        %% Make ind NrofCells based on header and read 36 values into each.
        Qstore = reshape(complex(R,I),[36,NrOfSarCells]);
        clear R I data header fname1 fid;
        Qavg = zeros(NrofCells,NrofCoils,NrofCoils);
         
        parfor k=1:NrofCells
        
           Qtm = squeeze(Qstore(:,k));
                       
            %% Put them as upper triangular matrices.
            Qtmfu = triu(ones(NrofCoils),0);
            Qtmfu(Qtmfu==1) = Qtm;
            Qavg(k,:,:) = (triu(Qtmfu,1))' + Qtmfu; %adding upper and lower matrices to form the full matrix, Hermetian.
       
       end
   
          %% Create output structure with details.
           Q.avg = Qavg;
           Q.index= index;
           Q.dimx = dimx;
           Q.dimy = dimy;
           Q.dimz = dimz;
           Q. NrOfSarCells =  NrOfSarCells;
           
            disp('Done reading index files and creating Q structure....');
           
           % Clean up the function local variables anyway
           clear Qavg index Qstore dimx dimy dimz  NrOfSarCells data fid fname NrofCoils NrofCells;
        
        
        
        
end
    