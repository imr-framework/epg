function Epwr = gen_E12ptQ(Ex,Ey,Ez,X,SigmabyRhox)
SigmabyRhoy=SigmabyRhox;
SigmabyRhoz=SigmabyRhox;
dbstop if error;
%% Get co-ordinates
if(numel(X)==1)
        M = size(Ex,1);
        N = size(Ex,2);
        P = size(Ex,3);
        [x,y,z]=ind2sub([M,N,P],X);
        X = [x;y;z];
else
        x=X(1);
        y=X(2);
        z=X(3);
end

%% Define co-ordinates of the 12-point cube formulation
X1 = [x,y+1,z];
X2 = [x,y,z+1];
X3 = [x,y+1,z+1];

%%
Y1 = [x+1,y,z];
Y2 = [x,y,z+1];
Y3 = [x+1,y,z+1];

%%
Z1 = [x+1,y,z];
Z2 = [x,y+1,z];
Z3 = [x+1,y+1,z];

%%
Ex1 = get_E(squeeze(Ex(X(1),X(2),X(3),:)),squeeze(SigmabyRhox(X(1),X(2),X(3))));
Ey1 = get_E(squeeze(Ey(X(1),X(2),X(3),:)),squeeze(SigmabyRhoy(X(1),X(2),X(3))));
Ez1 = get_E(squeeze(Ez(X(1),X(2),X(3),:)),squeeze(SigmabyRhoz(X(1),X(2),X(3))));
 
Expwr = Ex1 + get_E(squeeze(Ex(X1(1),X1(2),X1(3),:)),squeeze(SigmabyRhox(X1(1),X1(2),X1(3)))) + get_E(squeeze(Ex(X2(1),X2(2),X2(3),:)),squeeze(SigmabyRhox(X2(1),X2(2),X2(3)))) + get_E(squeeze(Ex(X3(1),X3(2),X3(3),:)),squeeze(SigmabyRhox(X3(1),X3(2),X3(3))));
Eypwr = Ey1 + get_E(squeeze(Ey(Y1(1),Y1(2),Y1(3),:)),squeeze(SigmabyRhoy(Y1(1),Y1(2),Y1(3)))) + get_E(squeeze(Ey(Y2(1),Y2(2),Y2(3),:)),squeeze(SigmabyRhoy(Y2(1),Y2(2),Y2(3)))) + get_E(squeeze(Ey(Y3(1),Y3(2),Y3(3),:)),squeeze(SigmabyRhoy(Y3(1),Y3(2),Y3(3))));
Ezpwr = Ez1 + get_E(squeeze(Ez(Z1(1),Z1(2),Z1(3),:)),squeeze(SigmabyRhoz(Z1(1),Z1(2),Z1(3)))) + get_E(squeeze(Ez(Z2(1),Z2(2),Z2(3),:)),squeeze(SigmabyRhoz(Z2(1),Z2(2),Z2(3)))) + get_E(squeeze(Ez(Z3(1),Z3(2),Z3(3),:)),squeeze(SigmabyRhoz(Z3(1),Z3(2),Z3(3))));
Epwr =0.125.*(Expwr + Eypwr + Ezpwr); %actually 0.25 but 0.25/2*Rho not captured by the SigmabyRho matrix


function Epwr =get_E(E,SigmabyRho)
Epwr = SigmabyRho.*(E*E');
