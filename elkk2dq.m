function [el_k,el_f] = elkk2dq(ntype,npe,fc,nfunc,el_x,el_y,nnod)
 
% Gauss points array
%gauss = [0.0,-1/sqrt(3),-sqrt(3/5),-0.86113631;0.0,1/sqrt(3),0.0,...
   % -0.33998104;0.0,0.0,sqrt(3/5),0.33998104;0.0,0.0,0.0,0.86113631];
% Gauss weights array

%weight = [2.0,1.0,5/9,0.34785485;0.0,1.0,8/9,0.65214515;...
   % 0.0,0.0,5/9,0.65214515;0.0,0.0,0.0,0.34785485];
W=1.0;
gauss=[sqrt(1/3),-sqrt(1/3),-sqrt(1/3),sqrt(1/3);sqrt(1/3),sqrt(1/3),...
    -sqrt(1/3),-sqrt(1/3)];

if ntype==2
    ngp=4;
% npea gives the size of the of element matrix   
    npea = 4;
else
    ngp=npe;
    npea=npe;
end


el_f = zeros(npea,1);
el_k = zeros(npea,npea);

for ngi =1:ngp
    for ngj=1:ngp
    xi = gauss(ngi,ngp); 
    eta=gauss(ngj,ngp);
    
    [sf,dsfxi,dsfeta,rj(1,1),rj(1,2),rj(2,1),rj(2,2)]=shape2dq(nfunc,npea,npe,xi,eta,el_x,el_y,nnod);
    
    
    %ax = ac(ng,1)+ac(ng,2)*x+ac(ng,3)*x*x;
    %cx = cc(ng,1)+cc(ng,2)*x+cc(ng,3)*x*x;
    fx = fc(ng,1)+fc(ng,2)*x+fc(ng,3)*x*x;
    %bx = bc(ng,1)+bc(ng,2)*x+fc(ng,3)*x*x;
    ex = ec(ngi,1)+ec(ngi,2)+ec(ngi,3);
    ey = ec(ngj,1)+ec(ngj,2)+ec(ngj,3);
%
    for i=1:npea
        el_f(i) = el_f(i)+fx*sf(i)*weight(ng,ngp)*rj;
    end
    
    if ntype==1
        for i=1:npea
            for j=1:npea
                el_k(i,j) = el_k(i,j)+(ex*((1/rj(1,1)*dsfxi(i)+1/rj(1,2)*dsfeta(i))...
                    *(1/rj(1,1)*dsfxi(j)+1/rj(1,2)*dsfeta(j)))+...
                    ey*((1/rj(2,1)*dsfxi(i)+1/rj(2,2)*dsfeta(i))*(1/rj(2,1)...
                    *dsfxi(j)+1/rj(2,2)*dsfeta(j))));
               
            end
        end
    end
end
    

end


