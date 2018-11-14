function [Qavg_df] = get_Qavg(Mass_cell,Mdef,Q,ind,s,t,u)

%% Variable initialization
Mass_air = 1.625e-7 + 1e-9;
D= size(Mass_cell);
Max_range =10; %growing the size of the cube, can be even more restrictive

%% 
          

[M,N,P,~,~]=size(Q);   
Qavg = zeros(length(ind),8,8); %Remove hardcoding for coils later
AvgSize = zeros(length(ind),1);
          
%%          
parfor k1=1:length(ind)

    [x,y,z]=ind2sub([M,N,P],ind(k1));
%     
       
%        if((Mass_cell(x,y,z) > Mass_air) && (chk==1))
          if((Mass_cell(x,y,z) > Mass_air))
              chk = chkcubair_global(2,Mass_cell,Mass_air,x,y,z);
              if(chk==1)
                Mcore = Mass_cell(x,y,z);

            Mcorner =0;Qcorner=0;
            Medge=0;Qedge=0;
            Mface=0;Qface=0;
            dim =0;
            Mtotal = Mcore + Mcorner + Medge + Mface;
            
%             Qtotal = Qcore + Qcorner + Qedge + Qface;    
                while((Mtotal < Mdef) && (dim < Max_range))
                        dim = dim +1;
                        Mcore = Mtotal;
  
                        %% Corner Mass
                        Mcorner = Mass_cell(x-dim,y-dim,z-dim) + Mass_cell(x-dim,y-dim,z+dim) + Mass_cell(x-dim,y+dim,z-dim) + Mass_cell(x-dim,y+dim,z+dim)+ ...
                                  Mass_cell(x+dim,y-dim,z-dim)+ Mass_cell(x+dim,y-dim,z+dim) + Mass_cell(x+dim,y+dim,z-dim) + Mass_cell(x+dim,y+dim,z+dim); 
                        
                        Qcorner = (squeeze(Q(x-dim,y-dim,z-dim,:,:)).* Mass_cell(x-dim,y-dim,z-dim) + ...
                                   squeeze(Q(x-dim,y-dim,z+dim,:,:)).* Mass_cell(x-dim,y-dim,z+dim) + ...
                                   squeeze(Q(x-dim,y+dim,z-dim,:,:)).* Mass_cell(x-dim,y+dim,z-dim) + ...
                                   squeeze(Q(x-dim,y+dim,z+dim,:,:)).* Mass_cell(x-dim,y+dim,z+dim)+ ...
                                   squeeze(Q(x+dim,y-dim,z-dim,:,:)).* Mass_cell(x+dim,y-dim,z-dim)+ ... 
                                   squeeze(Q(x+dim,y-dim,z+dim,:,:)).* Mass_cell(x+dim,y-dim,z+dim)+ ...
                                   squeeze(Q(x+dim,y+dim,z-dim,:,:)).* Mass_cell(x+dim,y+dim,z-dim)+ ...
                                   squeeze(Q(x+dim,y+dim,z+dim,:,:)).* Mass_cell(x+dim,y+dim,z+dim));       
                              
                              
                              
                              
                        %% Edge Mass
                        Medge=0;Qedge=0;
                        for k=-dim+1:dim-1 %Will start only from 2 onwards - clever imp from Philips
                            Medge = Medge + ...
                                    Mass_cell(x-k,y-dim,z-dim)+ Mass_cell(x-k,y-dim,z+dim)+ Mass_cell(x-k,y+dim,z-dim) + Mass_cell(x-k,y+dim,z+dim) +...
                                    Mass_cell(x-dim,y-k,z-dim) + Mass_cell(x-dim,y-k,z+dim) + Mass_cell(x+dim,y-k,z-dim)+ Mass_cell(x+dim,y-k,z+dim) +...
                                    Mass_cell(x-dim,y-dim,z-k) + Mass_cell(x-dim,y+dim,z-k) + Mass_cell(x+dim,y-dim,z-k) + Mass_cell(x+dim,y+dim,z-k);
                            
                            Qedge = Qedge + ...
                                    Q(x-k,y-dim,z-dim,:,:).* Mass_cell(x-k,y-dim,z-dim)+ ...
                                    Q(x-k,y-dim,z+dim,:,:).* Mass_cell(x-k,y-dim,z+dim)+ ...
                                    Q(x-k,y+dim,z-dim,:,:).* Mass_cell(x-k,y+dim,z-dim)+ ...
                                    Q(x-k,y+dim,z+dim,:,:).* Mass_cell(x-k,y+dim,z+dim)+ ...
                                    Q(x-dim,y-k,z-dim,:,:).* Mass_cell(x-dim,y-k,z-dim)+ ...
                                    Q(x-dim,y-k,z+dim,:,:).* Mass_cell(x-dim,y-k,z+dim)+ ...
                                    Q(x+dim,y-k,z-dim,:,:).* Mass_cell(x+dim,y-k,z-dim)+ ...
                                    Q(x+dim,y-k,z+dim,:,:).* Mass_cell(x+dim,y-k,z+dim)+ ...
                                    Q(x-dim,y-dim,z-k,:,:).* Mass_cell(x-dim,y-dim,z-k)+ ...
                                    Q(x-dim,y+dim,z-k,:,:).* Mass_cell(x-dim,y+dim,z-k)+ ...
                                    Q(x+dim,y-dim,z-k,:,:).* Mass_cell(x+dim,y-dim,z-k)+ ...
                                    Q(x+dim,y+dim,z-k,:,:).* Mass_cell(x+dim,y+dim,z-k);
    
                                

                        end


                        %% Face Mass
                        Mface=0; Qface=0;
                        for a = -dim+1:dim-1
                            for b=-dim+1:dim-1
                                Mface = Mface + ...
                                        Mass_cell(x+a,y+b,z+dim) + Mass_cell(x+a,y+b,z-dim)+ ...
                                        Mass_cell(x+a,y+dim,z+b) + Mass_cell(x+a,y-dim,z-b) + ...
                                        Mass_cell(x+dim,y+a,z+b) + Mass_cell(x-dim,y-a,z-b);
                                    
                                    
                                Qface = Qface + ...
                                        Q(x+a,y+b,z+dim,:,:).* Mass_cell(x+a,y+b,z+dim) + ...
                                        Q(x+a,y+b,z-dim,:,:).* Mass_cell(x+a,y+b,z-dim)+ ...
                                        Q(x+a,y+dim,z+b,:,:).* Mass_cell(x+a,y+dim,z+b)+ ...
                                        Q(x+a,y-dim,z-b,:,:).* Mass_cell(x+a,y-dim,z-b)+ ...
                                        Q(x+dim,y+a,z+b,:,:).* Mass_cell(x+dim,y+a,z+b)+ ...
                                        Q(x-dim,y-a,z-b,:,:).* Mass_cell(x-dim,y-a,z-b);
                            end
                        end


                       Mtotal = Mcore + Mcorner + Medge + Mface;
%                        Q_prev = Qtotal;
%                        Qtotal = Qcore + Qcorner + Qedge + Qface;

                end
                
                %% Add interpolation part
                kdf = Mcore - Mdef;       
                f=roots([Mcorner,Medge,Mface,kdf]);
                f = f(f>0);
       
                %% Get Qcore from previous cube
                r = dim-1;
                Qcore = 0;
                for i = x-r:x+r
                    for j= y-r:y+r
                        for k = z-r:z+r
                   
                              Qcore = Qcore + (squeeze(Q(i,j,k,:,:)).*Mass_cell(x,y,z));
                            
                        end
                    end
                end
                
%                 [Qt]= getQs(Q,Mass_cell,dim,x,y,z,f,Mdef);

                if((f>1)||(f<0))
                    error('value above expected range');
                end
                Qavg(k1,:,:) = ((squeeze(Qcorner).*(f^3)) + (squeeze(Qedge).*(f^2)) + (squeeze(Qface).*f) + (squeeze(Qcore)))./Mdef;
                
%                 figure;imagesc(abs((Qavg)));
%                 clear Qcore Mcore Qcorner Mcorner Qedge Medge Qface Mface;
                AvgSize(k1)= dim;   
              end
        else
            AvgSize(k1)=-1;                
        end
        
        
end

              Qavg_df= zeros(D,8,8);
              for   k=1:length(ind)
              Qavg_df(s(k),t(k),u(k),:,:) = squeeze(Qavg(k,:,:));
              end



           

    

  