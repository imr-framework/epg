function [status] = write_qmat(Q,SAR_type)

switch SAR_type
    
    case 'Global'
         [filename, pathname] = uiputfile('globalQMatrix_imp.tri', 'Save global Q matrix as');
         fid = fopen(fullfile(pathname,filename),'w+','l');
         
         %% Write header information in uint16
         NrofCoils = size(Q.Qtmf,1); %8
         NrofSarCells = 3; %body, head and exposed.
         
         fwrite(fid,NrofCoils,'uint32');
         fwrite(fid,NrofCoils,'uint32');
         fwrite(fid,NrofSarCells,'uint32');
        
         
         %% Write Q matrices in the order of body, head and exposed
         write_uppertri(Q.Qtmf,fid);
         write_uppertri(Q.Qhmf,fid);
         write_uppertri(Q.Qemf,fid);
         status=0;
         fclose all;
        
         
         
    case 'Local'
         [df_filename, df_pathname] = uigetfile('*.df', 'Open Bodysegmentation.df');
         Tissue_types= df_read(fullfile(df_pathname,df_filename));
         
         
         [tri_filename, tri_pathname] = uiputfile('avgQMatrix_imp.tri', 'Save local Q matrix as');
         tri_fid = fopen(fullfile(tri_pathname,tri_filename),'w+','l');
         
         [index_filename, index_pathname] = uiputfile('avgQMatrix_imp.index', 'Save local Q matrix indices as');
         index_fid = fopen(fullfile(index_pathname,index_filename),'w+','l');
         
         %% Write header information in uint16
         Qavg = Q;clear Q;
         [M,N,P,~,NrofCoils]=size(Qavg);
         Qavg = reshape(Qavg,[M*N*P,NrofCoils,NrofCoils]);
         S = Qavg>0;
         ind = find(squeeze(S(:,4,4)));
         NrofSarCells = length(ind);
         
         
         fwrite(tri_fid,NrofCoils,'uint32');
         fwrite(tri_fid,NrofCoils,'uint32');
         fwrite(tri_fid,NrofSarCells,'uint32');
                  
         fwrite(index_fid,M,'uint32');
         fwrite(index_fid,N,'uint32');
         fwrite(index_fid,P,'uint32');
         fwrite(index_fid,NrofSarCells,'uint32');
                  
         
         
         
         for k=1:length(ind)
             [x,y,z]=ind2sub([M,N,P],ind(k));
             label = Tissue_types(x,y,z);
             write_uppertri(squeeze(Qavg(ind(k),:,:)),tri_fid);  
             fwrite(index_fid,x-1,'uint16'); %compensate for the indexing in C
             fwrite(index_fid,y-1,'uint16');
             fwrite(index_fid,z-1,'uint16');
             fwrite(index_fid,label,'uint16');
         end
         
         
        status=0;
end

    function  write_uppertri(Qin,fid)
         Q = triu(Qin);
         Q_tri = Q(abs(Q)>0);
         if(length(Q_tri)~=36)
             error('Missing values in this matrix for 8 coils');
         end
         for k=1:length(Q_tri)
              fwrite(fid,real(Q_tri(k)),'float');
              fwrite(fid,imag(Q_tri(k)),'float');
         end
         

