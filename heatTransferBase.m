function Qfinb=heatTransferBase(u,kappa,nodes,elements)
%Heat transfer at base (per meter depth) in W/m
IEN=elements';
ref=1:size(nodes,1);
idxb=(nodes(:,1)<eps);
nbase=ref(idxb);
elembase=[];
for e=1:size(elements,1)
    for j=1:4
        jp1=mod(j,4)+1;
        edg=nnz(abs(nbase-elements(e,j))<eps)...
            +nnz(abs(nbase-elements(e,jp1))<eps);
        if edg==2
            n1=elements(e,j);
            n2=elements(e,jp1);
            elembase=[elembase;[e,n1,n2]];
        end
    end
end
Qfinb=0;
for i=1:size(elembase,1)
    e=elembase(i,1); v1=elembase(i,2); v2=elembase(i,3);
    idx=[find(abs(IEN(:,e)-v1)<eps),find(abs(IEN(:,e)-v2)<eps)];
    xnod=nodes(IEN(:,e),:);
    aloc=u(IEN(:,e));
    nml=[-1;0];%normal vector
    n1D=2;
    [s,w1D]=mylegendrepts(n1D);
    switch idx(1)
    case 1
        xi2D=[s,-ones(n1D,1)];
    case 2
        xi2D=[ones(n1D,1),s];
    case 3
        xi2D=[-s,ones(n1D,1)];
    case 4
        xi2D=[-ones(n1D,1),-s];
    end
    for l=1:length(w1D)
        [phi,gradxiphi]=shapeQuad(xi2D(l,:));
        [~,~,gradxphi]=shapeQuadElem(xnod,phi,gradxiphi);
        gradxphiedgenml=nml'*gradxphi*aloc;
        detJedge=norm(xnod(idx(1),:)-xnod(idx(2),:))/2;
        Qfinb=Qfinb+kappa*gradxphiedgenml*detJedge*w1D(l);
    end
end