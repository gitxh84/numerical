function  fheat = HeatSupply(r,xnod)
% r: function handle of heat supply per unit volume
% xnod: 4 * 2 matrix of coordinates
% rtype:
%   fheat: 4 * 1

    fheat = zeros(4,1);
    
    n1D = 2;
    [xi2D,w2D] = GaussLeg2DQuad(n1D);  
    
    for a = 1:4 
        temp = 0;
        for i = 1:(n1D^2)
            xieta = xi2D(i,:);
            [phi, gradxiphi] = shapeQuad(xieta);
            [xy, detJ, gradxphi] = shapeQuadElem(xnod, phi, gradxiphi);

            temp = temp + w2D(1,i) * (r(xy(1),xy(2)) * phi(1,a) * detJ);
        end
        
        fheat(a,1) = temp;
                
    end

end

