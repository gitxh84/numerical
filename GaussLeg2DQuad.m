function [xi2D,w2D] = GaussLeg2DQuad(n1D)
% n1D: number of 1D Gauss-Legendre points in x and y
% rtype: 
%   xi2D: n1D^2 * 2
%   w2D: 1 * n1D^2

    n_int = n1D ^ 2;
    xi2D = zeros(n_int, 2);
    w2D = zeros(1, n_int);
    
    [x,w_x] = mylegendrepts(n1D);
    [y,w_y] = mylegendrepts(n1D);

    for l = 0:(n_int - 1)
        
        i = floor(l / n1D) + 1;
        j = mod(l, n1D) + 1;
        
        xi2D(l+1,1) = x(i,1);
        xi2D(l+1,2) = y(j,1);
        w2D(1, l+1) = w_x(1,i) * w_y(1,j);

    end

end

