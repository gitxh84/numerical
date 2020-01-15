clear all
%Geometry
tb=0.001;%in m
L=0.005;%in m
h=tb/32;%rough element size in m
tf=0.0004;%in m
[elementsR,nodesR]=rectangularFinMeshGen(tb,L,h);
[elementsC,nodesC]=choppedTriFinMeshGen(tb,tf,L,h);
[elementsP,nodesP]=parabolicTriFinMeshGen(tb,tf,L,h);
figure
subplot(3,1,1);plot(nodesR(:,1),nodesR(:,2),'o');axis equal;axis tight;
subplot(3,1,2);plot(nodesC(:,1),nodesC(:,2),'o');axis equal;axis tight;
subplot(3,1,3);plot(nodesP(:,1),nodesP(:,2),'o');axis equal;axis tight;
%Boundary conditions
Tb=350; %temperature at base in K
Troom=300; %room temperature in K
kappa=205; %Aluminum thermal conductivity in W/(mK)
hconv=25; %Heat transfer coefficient at tip of fin
r=@(x,y) 0;








% SCENARIO 1: DIRICHLET BOUNDARY CONDITIONS
%First lets check the fins with full temperature (Dirichlet) BCs
[dirichlet]=FinFullDirichletBCs(nodesR,L,tb,Tb,Troom,kappa,hconv);
% FIN 1: Rectangular Fin
uhR=FE2DNonzeroDirichlet(nodesR,elementsR,dirichlet,kappa,r);
uavgR=Taverage(uhR,nodesR,elementsR);
%Heat transfer at base (per meter depth) in W/m
QfinbR=heatTransferBase(uhR,kappa,nodesR,elementsR);
%Efficiency and Effectiveness of fin
LfinR=2*L+tb;
FinEffecR=QfinbR/(hconv*tb*(Tb-Troom));
FinEfficR=FinEffecR*(tb/LfinR);
% FIN 2: Rectangular Fin
uhC=FE2DNonzeroDirichlet(nodesC,elementsC,dirichlet,kappa,r);
uavgC=Taverage(uhC,nodesC,elementsC);
%Heat transfer at base (per meter depth) in W/m
QfinbC=heatTransferBase(uhC,kappa,nodesC,elementsC);
%Efficiency and Effectiveness of fin
LfinC=2*sqrt(((tb-tf)/2)^2+L^2)+tf;
FinEffecC=QfinbC/(hconv*tb*(Tb-Troom));
FinEfficC=FinEffecC*(tb/LfinC);
% FIN 3: Parabolic Fin
uhP=FE2DNonzeroDirichlet(nodesP,elementsP,dirichlet,kappa,r);
uavgP=Taverage(uhP,nodesP,elementsP);
%Heat transfer at base (per meter depth) in W/m
QfinbP=heatTransferBase(uhP,kappa,nodesP,elementsP);
%Efficiency and Effectiveness of fin
LfinP=2*(0.5*sqrt((tf-tb)^2+L^2)+(L^2/(2*(tf-tb)))*asinh((tf-tb)/L))+tf;
FinEffecP=QfinbP/(hconv*tb*(Tb-Troom));
FinEfficP=FinEffecP*(tb/LfinP);
% SUMMARY OF RESULTS SCENARIO 1: each column is a fin
FinsPerformDir=[uavgR,uavgC,uavgP;...
                   FinEffecR,FinEffecC,FinEffecP;...
                   FinEfficR,FinEfficC,FinEfficP];

               
               

               
               
               
               
  
% SCENARIO 2: DIRICHLET+NEUMANN BOUNDARY CONDITIONS
%Lets check the fins w temperature (Dirichlet) BCs and heat (Neumann) BCs
[dirichlet,neumann]=FinMixedNeumannBCs(nodesR,elementsR,L,tb,Tb,Troom,kappa,hconv);
% FIN 1: Rectangular Fin
uhR=FE2DDirNeu(nodesR,elementsR,dirichlet,neumann,kappa,r);
uavgR=Taverage(uhR,nodesR,elementsR);
%Heat transfer at base (per meter depth) in W/m
QfinbR=heatTransferBase(uhR,kappa,nodesR,elementsR);
%Efficiency and Effectiveness of fin
LfinR=2*L+tb;
FinEffecR=QfinbR/(hconv*tb*(Tb-Troom));
FinEfficR=FinEffecR*(tb/LfinR);
% FIN 2: Rectangular Fin
uhC=FE2DDirNeu(nodesC,elementsC,dirichlet,neumann,kappa,r);
uavgC=Taverage(uhC,nodesC,elementsC);
%Heat transfer at base (per meter depth) in W/m
QfinbC=heatTransferBase(uhC,kappa,nodesC,elementsC);
%Efficiency and Effectiveness of fin
LfinC=2*sqrt(((tb-tf)/2)^2+L^2)+tf;
FinEffecC=QfinbC/(hconv*tb*(Tb-Troom));
FinEfficC=FinEffecC*(tb/LfinC);
% FIN 3: Parabolic Fin
uhP=FE2DDirNeu(nodesP,elementsP,dirichlet,neumann,kappa,r);
uavgP=Taverage(uhP,nodesP,elementsP);
%Heat transfer at base (per meter depth) in W/m
QfinbP=heatTransferBase(uhP,kappa,nodesP,elementsP);
%Efficiency and Effectiveness of fin
LfinP=2*(0.5*sqrt((tf-tb)^2+L^2)+(L^2/(2*(tf-tb)))*asinh((tf-tb)/L))+tf;
FinEffecP=QfinbP/(hconv*tb*(Tb-Troom));
FinEfficP=FinEffecP*(tb/LfinP);
% SUMMARY OF RESULTS SCENARIO 2: each column is a fin
FinsPerformDirNeu=[uavgR,uavgC,uavgP;...
                   FinEffecR,FinEffecC,FinEffecP;...
                   FinEfficR,FinEfficC,FinEfficP];

               
               
               
               
               
               
               
               
               
           
               
               
               
               
% SCENARIO 3: DIRICHLET+ROBIN BOUNDARY CONDITIONS: most realistic
%Lets check the fins w temperature (Dirichlet) BCs and convection (Robin) BCs
[dirichlet,neumann,robin]=FinMixedRobinBCs(nodesR,elementsR,L,tb,Tb,Troom,kappa,hconv);
% FIN 1: Rectangular Fin
uhR=FE2DDirNeuRob(nodesR,elementsR,dirichlet,neumann,robin,kappa,hconv,r);
uavgR=Taverage(uhR,nodesR,elementsR);
%Heat transfer at base (per meter depth) in W/m
QfinbR=heatTransferBase(uhR,kappa,nodesR,elementsR);
%Efficiency and Effectiveness of fin
LfinR=2*L+tb;
FinEffecR=QfinbR/(hconv*tb*(Tb-Troom));
FinEfficR=FinEffecR*(tb/LfinR);
% FIN 2: Rectangular Fin
uhC=FE2DDirNeuRob(nodesC,elementsC,dirichlet,neumann,robin,kappa,hconv,r);
uavgC=Taverage(uhC,nodesC,elementsC);
%Heat transfer at base (per meter depth) in W/m
QfinbC=heatTransferBase(uhC,kappa,nodesC,elementsC);
%Efficiency and Effectiveness of fin
LfinC=2*sqrt(((tb-tf)/2)^2+L^2)+tf;
FinEffecC=QfinbC/(hconv*tb*(Tb-Troom));
FinEfficC=FinEffecC*(tb/LfinC);
% FIN 3: Parabolic Fin
uhP=FE2DDirNeuRob(nodesP,elementsP,dirichlet,neumann,robin,kappa,hconv,r);
uavgP=Taverage(uhP,nodesP,elementsP);
%Heat transfer at base (per meter depth) in W/m
QfinbP=heatTransferBase(uhP,kappa,nodesP,elementsP);
%Efficiency and Effectiveness of fin
LfinP=2*(0.5*sqrt((tf-tb)^2+L^2)+(L^2/(2*(tf-tb)))*asinh((tf-tb)/L))+tf;
FinEffecP=QfinbP/(hconv*tb*(Tb-Troom));
FinEfficP=FinEffecP*(tb/LfinP);
% SUMMARY OF RESULTS SCENARIO 3: each column is a fin
FinsPerformDirRob=[uavgR,uavgC,uavgP;...
                   FinEffecR,FinEffecC,FinEffecP;...
                   FinEfficR,FinEfficC,FinEfficP];


               
            
               
               
               
               
               
               
           
               
               
               
               
               
               
               
function [dirichlet]=FinFullDirichletBCs(nodesR,L,tb,Tb,Troom,kappa,hconv)
ref=1:size(nodesR,1);
idxDir=(nodesR(:,1)<eps|nodesR(:,1)>L-eps|nodesR(:,2)<-tb/2+eps|nodesR(:,2)>tb/2-eps);
nDir=ref(idxDir);
%We choose BCs independent of y as they will be good for all fin
%configurations. It is possible to choose more complex BCs.
m=sqrt(2*hconv/(kappa*tb));
gD=@(x,y) Troom+((Tb-Troom)/cosh(m*L))*cosh(m*(L-x)); %Temperature
dirichlet=[nDir',gD(nodesR(nDir',1),nodesR(nDir',2))];
end

function [dirichlet,neumann]=FinMixedNeumannBCs(nodesR,elements,L,tb,Tb,Troom,kappa,hconv)
ref=1:size(nodesR,1);
idxDir=(nodesR(:,1)<eps);
idxNeu=(nodesR(:,1)>L-eps|nodesR(:,2)<-tb/2+eps|nodesR(:,2)>tb/2-eps);
nDir=ref(idxDir);
nNeu=ref(idxNeu);
%We choose BCs independent of y as they will be good for all fin
%configurations. It is possible to choose more complex BCs.
m=sqrt(2*hconv/(kappa*tb));
gD=@(x,y) 0*x+Tb; %Temperature
dirichlet=[nDir',gD(nodesR(nDir',1),nodesR(nDir',2))];
gN=@(x,y) -hconv*((Tb-Troom)/cosh(m*L))*cosh(m*(L-x));
neumann=[];
for e=1:size(elements,1)
    for j=1:4
        jp1=mod(j,4)+1;
        edg=nnz(abs(nNeu-elements(e,j))<eps)...
            +nnz(abs(nNeu-elements(e,jp1))<eps);
        if edg==2
            n1=elements(e,j);
            n2=elements(e,jp1);
            gN1=gN(nodesR(n1,1),nodesR(n1,2));
            gN2=gN(nodesR(n2,1),nodesR(n2,2));
            neumann=[neumann;[e,n1,n2,gN1,gN2]];
        end
    end
end
end

function [dirichlet,neumann,robin]=FinMixedRobinBCs(nodesR,elements,L,tb,Tb,Troom,kappa,hconv)
ref=1:size(nodesR,1);
idxDir=(nodesR(:,1)<eps);
idxRob=(nodesR(:,1)>L-eps|nodesR(:,2)<-tb/2+eps|nodesR(:,2)>tb/2-eps);
nDir=ref(idxDir);
nRob=ref(idxRob);
%We choose BCs independent of y as they will be good for all fin
%configurations. It is possible to choose more complex BCs.
gD=@(x,y) 0*x+Tb; %Temperature
dirichlet=[nDir',gD(nodesR(nDir',1),nodesR(nDir',2))];
neumann=[];
gR=@(x,y) Troom;
robin=[];
for e=1:size(elements,1)
    for j=1:4
        jp1=mod(j,4)+1;
        edg=nnz(abs(nRob-elements(e,j))<eps)...
            +nnz(abs(nRob-elements(e,jp1))<eps);
        if edg==2
            n1=elements(e,j);
            n2=elements(e,jp1);
            gR1=gR(nodesR(n1,1),nodesR(n1,2));
            gR2=gR(nodesR(n2,1),nodesR(n2,2));
            robin=[robin;[e,n1,n2,gR1,gR2]];
        end
    end
end
end

function [elements,nodes]=parabolicTriFinMeshGen(tb,tf,L,h)
[elements,nodes]=rectangularFinMeshGen(tb,L,h);
ctrans=@(x,y) [x,(2/tb)*(((tf-tb)/(2*L^2))*x.^2+tb/2).*y];
nodes=ctrans(nodes(:,1),nodes(:,2));
end

function [elements,nodes]=choppedTriFinMeshGen(tb,tf,L,h)
[elements,nodes]=rectangularFinMeshGen(tb,L,h);
xnod=[0,-tb/2;L,-tf/2;L,tf/2;0,tb/2]; %Nodes of fin
phi1=@(x,y) 0.25*(1-(2/L)*(x-(L/2))).*(1-2*y/tb);
phi2=@(x,y) 0.25*(1+(2/L)*(x-(L/2))).*(1-2*y/tb);
phi3=@(x,y) 0.25*(1+(2/L)*(x-(L/2))).*(1+2*y/tb);
phi4=@(x,y) 0.25*(1-(2/L)*(x-(L/2))).*(1+2*y/tb);
trans=@(x,y) phi1(x,y)*xnod(1,:)+phi2(x,y)*xnod(2,:)+...
                phi3(x,y)*xnod(3,:)+phi4(x,y)*xnod(4,:);
nodes=trans(nodes(:,1),nodes(:,2));
end

function [elements,nodes]=rectangularFinMeshGen(tb,L,h)
Nx=ceil(L/h); hx=L/Nx;
Ny=ceil(tb/h); hy=tb/Ny;
[Y,X]=meshgrid(-tb/2:hy:tb/2,0:hx:L);
nodes=[X(:),Y(:)];
Nnodes=(Nx+1)*(Ny+1);
ref=1:Nnodes;
idx0=logical(mod(ref',(Nx+1)));
idx1=idx0;
idx1(Nnodes-Nx:Nnodes)=0;
idx2=circshift(idx0,1);
idx2(Nnodes-Nx:Nnodes)=0;
idx3=circshift(idx0,1);
idx3(1:Nx+1)=0;
idx4=idx0;
idx4(1:Nx+1)=0;
elements=zeros(Nx*Ny,4);
elements(:,1)=ref(idx1);
elements(:,2)=ref(idx2);
elements(:,3)=ref(idx3);
elements(:,4)=ref(idx4);
end

% 
% %Plotting
% trisurf(elements, nodes(:,1), nodes(:,2), uh, ...
% 'EdgeColor', 'interp', 'FaceColor', 'interp' );
% xlabel('x') ; ylabel('y');
% title('Temperature Plot (K)')
% colorbar('location','eastoutside')
