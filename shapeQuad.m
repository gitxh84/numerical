function [phi,gradxiphi] = shapeQuad(xieta)
% xieta: 1 * 2 coordinate vector
% rtype:
%   phi: 1 * 4
%   gradxiphi: 2 * 4

    xi = xieta(1,1);
    eta = xieta(1,2);
    
    phi = [(1/4)*(1-xi)*(1-eta),(1/4)*(1+xi)*(1-eta),(1/4)*(1+xi)*(1+eta),(1/4)*(1-xi)*(1+eta)];
    
    gradxiphi = zeros(2,4);
    gradxiphi(1,1) = eta/4 - 1/4;
    gradxiphi(1,2) = 1/4 - eta/4;
    gradxiphi(1,3) = eta/4 + 1/4;
    gradxiphi(1,4) = -eta/4 - 1/4;
    
    gradxiphi(2,1) = xi/4 - 1/4;
    gradxiphi(2,2) = -xi/4 - 1/4;
    gradxiphi(2,3) = xi/4 + 1/4;
    gradxiphi(2,4) = 1/4 - xi/4;
    
end

