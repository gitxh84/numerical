function L2err = L2error(uh,uex,nodes,elements)
% rtype: number of L2 error between uh and uex

    n1D = 2;
    n_int = n1D ^ 2;
    [xi2D,w2D] = GaussLeg2DQuad(n1D);
    n_el = size(elements,1);

    sum_of_square = 0;
    for e=1:n_el
        xnod = nodes(elements(e,:),:);
        alpha = uh(elements(e,:));
        
        temp = 0;
        for l = 1:n_int
            
            xieta = xi2D(l,:);
            [phi,gradxiphi]=shapeQuad(xieta);
            [xy,detJ,gradxphi]=shapeQuadElem(xnod,phi,gradxiphi);
            
            uex_stuff = uex(xy(1,1),xy(1,2));
            uh_stuff = 0;
            for p = 1:4
                uh_stuff = uh_stuff + alpha(p,1) * phi(1,p);
            end
            
            temp = temp + (uex_stuff - uh_stuff) ^ 2 * detJ * w2D(l);
            
        end
        
        sum_of_square = sum_of_square + temp;
        
    end

    L2err = sqrt(sum_of_square);


end

