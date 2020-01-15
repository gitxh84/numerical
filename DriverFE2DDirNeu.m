%Driver to test codes
L=1;
N=128;
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
idxDir=(nodes(:,1)>L-eps|nodes(:,2)<eps|nodes(:,2)>L-eps);
nDir=ref(idxDir);
gD=@(x,y) 0*x.*y;
dirichlet=[nDir',gD(nodes(nDir',1),nodes(nDir',2))];
idxNeu=(nodes(:,1)<eps);
nNeu=ref(idxNeu);
gN=@(x,y) -pi*cos(pi*x).*sin(pi*y);
neumann=[];
for e=1:size(elements,1)
    for j=1:4
        jp1=mod(j,4)+1;
        edg=nnz(abs(nNeu-elements(e,j))<eps)...
            +nnz(abs(nNeu-elements(e,jp1))<eps);
        if edg==2
            n1=elements(e,j);
            n2=elements(e,jp1);
            gN1=gN(nodes(n1,1),nodes(n1,2));
            gN2=gN(nodes(n2,1),nodes(n2,2));
            neumann=[neumann;[e,n1,n2,gN1,gN2]];
        end
    end
end
plot(nodes(:,1),nodes(:,2),'o')

kappa=1;
r=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
uh=FE2DDirNeu(nodes,elements,dirichlet,neumann,kappa,r);
uex=@(x,y) sin(pi*x).*sin(pi*y);
uexact=uex(nodes(:,1),nodes(:,2));

%Plotting
trisurf(elements, nodes(:,1), nodes(:,2), uh, ...
'EdgeColor', 'interp', 'FaceColor', 'interp' );
xlabel('x') ; ylabel('y');
title('Temperature Plot (K)')
colorbar('location','eastoutside')

err=(1/length(uh))*norm(uh-uexact);
L2err=L2error(uh,uex,nodes,elements);