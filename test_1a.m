n1D = 4;
n = 6;
m = 6;

[xi2D,w2D] = GaussLeg2DQuad(n1D);

f = @(xi,eta) (xi.^n).*(eta.^m);
numInt = w2D*f(xi2D(:,1),xi2D(:,2));
exactInt=(1/(n+1))*(1^(n+1)-(-1)^(n+1))*(1/(m+1))*(1^(m+1) - (-1)^(m+1));

disp(numInt);                               % 0.0816
disp(exactInt);                             % 0.0816
disp(abs(exactInt - numInt) < eps(1));      % 1 (True)