function uavg = Taverage(uh,nodes,elements)
% rtype: number of average temperature

    n1D = 2;
    n_int = n1D ^ 2;
    n_el = size(elements,1);
    [xi2D,w2D] = GaussLeg2DQuad(n1D);
    
    integral_of_1 = 0;
    integral_of_uh = 0;
    for e=1:n_el
        
        xnod = nodes(elements(e,:),:);
        alpha = uh(elements(e,:));
        
        temp_integral_of_1 = 0;
        temp_integral_of_uh = 0;
        for l = 1:n_int
            
            xieta = xi2D(l,:);
            [phi,gradxiphi]=shapeQuad(xieta);
            [xy,detJ,gradxphi]=shapeQuadElem(xnod,phi,gradxiphi);
            
            temp = 0;
            for a = 1:4
                temp = temp + alpha(a,1) * phi(1,a);
            end
            
            temp_integral_of_1 = temp_integral_of_1 + 1 * detJ * w2D(l);
            temp_integral_of_uh = temp_integral_of_uh + temp * detJ * w2D(l);

        end
        
        integral_of_1 = integral_of_1 + temp_integral_of_1;
        integral_of_uh = integral_of_uh + temp_integral_of_uh;
        
    end

    uavg = 1/integral_of_1 * integral_of_uh;


end

