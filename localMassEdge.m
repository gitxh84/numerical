function medge=localMassEdge(xnod)
n1D=2;
[s,w1D]=mylegendrepts(n1D);
medge=zeros(size(xnod,1));
for l=1:length(w1D)
    psiedge=[0.5*(1-s(l)),0.5*(1+s(l))];
    detJ=norm(xnod(1,:)-xnod(2,:))/2;
    medge=medge+psiedge'*psiedge*detJ*w1D(l);
end