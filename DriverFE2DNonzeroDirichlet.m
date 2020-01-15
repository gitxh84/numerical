%Driver to test codes
L=1;
N=8;
h=L/N;
[Y,X]=meshgrid(0:h:L);
nodes=[X(:),Y(:)];
ref=1:(N+1)^2;
idx0=logical(mod(ref',(N+1)));
idx1=idx0;
idx1((N+1)^2-N:(N+1)^2)=0;
idx2=circshift(idx0,1);
idx2((N+1)^2-N:(N+1)^2)=0;
idx3=circshift(idx0,1);
idx3(1:N+1)=0;
idx4=idx0;
idx4(1:N+1)=0;
elements=zeros(N^2,4);
elements(:,1)=ref(idx1);
elements(:,2)=ref(idx2);
elements(:,3)=ref(idx3);
elements(:,4)=ref(idx4);
idxDir=(nodes(:,1)<eps|nodes(:,1)>L-eps|nodes(:,2)<eps|nodes(:,2)>L-eps);
nDir=ref(idxDir);
gD=@(x,y) x.*y;
dirichlet=[nDir',gD(nodes(nDir',1),nodes(nDir',2))];
plot(nodes(:,1),nodes(:,2),'o')

kappa=1;
r=@(x,y) 0;
uh=FE2DNonzeroDirichlet(nodes,elements,dirichlet,kappa,r);
uex=@(x,y) x.*y;
uexact=uex(nodes(:,1),nodes(:,2));

%Plotting
trisurf(elements, nodes(:,1), nodes(:,2), uh, ...
'EdgeColor', 'interp', 'FaceColor', 'interp' );
xlabel('x') ; ylabel('y');
title('Temperature Plot (K)')
colorbar('location','eastoutside')

err=(1/length(uh))*norm(uh-uexact);
L2err=L2error(uh,uex,nodes,elements);