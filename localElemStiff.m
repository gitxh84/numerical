function ke = localElemStiff(xnod)
% xnod: 4 * 2 matrix of coordinates
% rtype:
%   ke: 4 * 4 local stiffness matrix
    
    ke = zeros(4,4);
    [xi2D,w2D] = GaussLeg2DQuad(2);
    
    xieta = xi2D(1,:);
    [phi,gradxiphi]=shapeQuad(xieta);
    [xy,detJ_1,gradxphi_1]=shapeQuadElem(xnod,phi,gradxiphi);
    
    xieta = xi2D(2,:);
    [phi,gradxiphi]=shapeQuad(xieta);
    [xy,detJ_2,gradxphi_2]=shapeQuadElem(xnod,phi,gradxiphi);
    
    xieta = xi2D(3,:);
    [phi,gradxiphi]=shapeQuad(xieta);
    [xy,detJ_3,gradxphi_3]=shapeQuadElem(xnod,phi,gradxiphi);
    
    xieta = xi2D(4,:);
    [phi,gradxiphi]=shapeQuad(xieta);
    [xy,detJ_4,gradxphi_4]=shapeQuadElem(xnod,phi,gradxiphi);
    
    for i =1:4
        for j = 1:4
            sum = 0;
           
            for m = 1:4
                if m == 1
                    detJ = detJ_1;
                    gradxphi = gradxphi_1;
                    
                elseif m == 2
                    detJ = detJ_2;
                    gradxphi = gradxphi_2;
                    
                elseif m == 3
                    detJ = detJ_3;
                    gradxphi = gradxphi_3;
                    
                elseif m == 4
                    detJ = detJ_4;
                    gradxphi = gradxphi_4;
                    
                end
               
                product = transpose(gradxphi(:,j))*gradxphi(:,i)*detJ*w2D(m);
                sum = sum + product;
                
            end
            
            ke(i,j) = sum;
            
        end
    end


end

