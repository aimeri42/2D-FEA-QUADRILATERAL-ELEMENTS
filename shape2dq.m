function [sf,dsfxi,dsfeta,rj(1,1),rj(1,2),rj(2,1),rj(2,2)] = shape2dq(nfunc,npea,npe,xi,eta,el_x,el_y,nnod)

sf = zeros(npe,1);
dsf = zeros(npe,1);
ddsf = zeros(npe,1);
gdsf = zeros(npe,1);
gddsf = zeros(npe,1);
%
if nfunc==1
    sf(1) = 1/4*(1-eta)*(1-xi);
    sf(2) = 1/4*(1-eta)*(1+xi);
    sf(3) = 1/4*(1+eta)*(1+xi);
    sf(4) = 1/4*(1+eta)*(1-xi);
    dsfxi(1)= -1/4*(1-eta);
    dsfxi(2)=1/4*(1-eta);
    dsfxi(3)=1/4*(1+eta);
    dsfxi(4)=-1/4*(1+eta);
    dsfeta(1)=-1/4*(1-xi);
    dsfeta(2)=-1/4*(1+xi);
    dsfeta(3)=1/4*(1+xi);
    dsfeta(4)=1/4*(1-xi);
    
    
else
    sf(1) = -(1/2)*xi*(1-xi);
    dsf(1) = (-1/2)+xi;
    sf(2) = 1-xi^2;
    dsf(2) = -2*xi;
    sf(3) = (1/2)*xi*(1+xi);
    dsf(3) = (1/2)+xi;
end

%    Jacobian

for i= 1:nnod
rj(1,1)=rj(1,1)+dsfxi(i)*el_x(i);
rj(1,2)=rj(1,2)+dsfxi(i)*el_y(i);
rj(2,1)=rj(2,1)+dsfeta(i)*el_x(i);
rj(2,2)=rj(2,2)+dsfeta(i)*el_y(i);

end

%	Shape functions derivatives in global coordinates
for i = 1:npea
    gdsf(i) = dsf(i)/rj;
end

end

