xieta = [0,0];
xnod = [-1,-1 ; 
        1,-2 ;
        2,2 ;
        -1,1];
    
[phi,gradxiphi] = shapeQuad(xieta);
[xy,detJ,gradxphi] = shapeQuadElem(xnod,phi,gradxiphi);

disp(xy);           % correct 
disp(detJ);         % correct
disp(gradxphi);     % correct