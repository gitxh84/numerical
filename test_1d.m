r = @(x,y) 1+x+y;
xnod = [-1,-1; 1,-1; 2,1; 0,1];
fheat = HeatSupply(r,xnod);

disp(fheat);
    