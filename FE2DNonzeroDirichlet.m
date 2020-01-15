function uh=FE2DNonzeroDirichlet(nodes,elements,dirichlet,kappa,r)
[ID,IEN,LM]=locator(nodes,elements,dirichlet);
n_el=size(IEN,2); %Number of elements (also n_el=size(elements,1))
K=zeros(max(ID));
F=zeros(max(ID),1);
%When LM(a,e)=0, Dir is convenient since the Dirichlet BC will be given by
%Dir(IEN(a,e)). It might be useful in implementing nonzero Dirichlet BCs.
Dir=sparse(dirichlet(:,1),1,dirichlet(:,2),size(nodes,1),1);
for e=1:n_el
    xnod=nodes(IEN(:,e),:);
    %...construct K and F with Dirichlet BCs
    
    ke = localElemStiff(xnod);
    fe = HeatSupply(r,xnod);
    for p = 1:4
        i = LM(p,e);
        if i > 0
            F(i) = F(i) + fe(p);
            for b = 1:4
                j = LM(b,e);
                if j > 0
                    K(i,j) = K(i,j) + kappa * ke(p,b);
                elseif j == 0
                    g = Dir(IEN(b,e));
                    F(i) = F(i) - g * kappa * ke(p,b);
                end
            end
        end
    end
    
    
    
    
    
    
end
alpha=K\F;
uh=zeros(size(nodes,1),1);
uh(~logical(ID))=dirichlet(:,2);
uh(logical(ID))=alpha;
    