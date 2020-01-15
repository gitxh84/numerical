function uh=FE2DZeroDirichlet(nodes,elements,dirichlet,kappa,r)
[ID,IEN,LM]=locator(nodes,elements,dirichlet);
n_el=size(IEN,2); %Number of elements (also n_el=size(elements,1))
K=zeros(max(ID));
F=zeros(max(ID),1);
for e=1:n_el
    xnod=nodes(IEN(:,e),:);
    %...construct K and F without any Dirichlet BCs
    
    fe = HeatSupply(r,xnod);
    ke = localElemStiff(xnod);   
    [num_row, num_col] = size(LM);
    
    for ind_1 = 1:num_row
        i = LM(ind_1,e);
        if i > 0
            F(i) = F(i) + fe(ind_1);

            for ind_2 = 1:num_row
                j = LM(ind_2,e);
                if j > 0
                   K(i,j) = K(i,j) + kappa * ke(ind_1,ind_2);
                end
            end

        end
    end

end
alpha=K\F;
uh=zeros(size(nodes,1),1);
uh(~logical(ID))=dirichlet(:,2);
uh(logical(ID))=alpha;
    