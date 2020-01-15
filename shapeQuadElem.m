function [xy,detJ,gradxphi] = shapeQuadElem(xnod,phi,gradxiphi)
% xnod: 4 * 2 matrix of coordinates
% phi: 1 * 4 output of shapeQuad
% gradxiphi: 2 * 4 output of shapeQuad
% rtype:
%   xy: 1 * 2
%   detJ: number
%   gradxphi: 2 * 4

    xy = zeros(1,2);
    for i = 1:4 
        xy = xy + phi(1,i) * xnod(i,:);
    end
   
    J = zeros(2,2);
    for i = 1:4
        J = J + gradxiphi(:,i) * xnod(i,:);
    end
    
    detJ = det(J);
    gradxphi = J\gradxiphi;
 
end

