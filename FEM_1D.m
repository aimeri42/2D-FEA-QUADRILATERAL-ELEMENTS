% ===================
% 1-D FEM
% ME6350: FINITE ELEMENT ANALYSIS
% Astrit Imeri - T00249444
% % ===================
% ncon - connectivity array
% ndfn - number of unknowns (d.o.f) per node
% nelem - number of elements 
% nnod - number of nodes
% nnbc, nebc: number of known NBC and EBC
% inbc,iebc: location of known NBC and EBC
% vnbc,vebc: values of known NBC and EBC
%ntype - 1: 2nd order ODE, 2: 4th order ODE
%nfunc - 1: linear 2: quadratic
%gk - global coefficient matrix, gf- global right hand side vector

clc;
% Call the input file
input_file


npe=2;
if ntype==2
    ndfn=3;
else
    ndfn=1;
    if nfunc==2
        npe=3;
    end
end

% End reading data

% Calculate number of total equations to solve

neq = nnod;
g_k = zeros(neq,neq);
g_f = zeros(neq,1);
el_x = zeros(npe,1);
el_y= zeros(npe,1);
% Start loop over number of elements to calculate the element matrices
% and assemble them into the global matrices
%

for n=1:nelem
    for i=1:npe
        el_x(i) = x(ncon(n,i));
        el_y(i)= y(ncon(n,i));
    end
    el_x;
    el_y;
    
    [el_k,el_f] = elkk2dq(ntype,npe,fc,nfunc,el_x,el_y,nnod);
    
    % Assembly
    for i = 1:npe
        nr = (ncon(n,i)-1)*ndfn;
        for ii = 1:ndfn
            nr = nr + 1;
            l = (i-1)*ndfn+ii;
            g_f(nr) = g_f(nr)+el_f(l);
            for j = 1:npe
                nc = (ncon(n,j)-1)*ndfn;
                for jj = 1:ndfn
                    m = (j-1)*ndfn+jj;
                    nc = nc + 1;
                    g_k(nr,nc) =g_k(nr,nc)+el_k(l,m); 
                end
            end
        end
    end
end

% Apply known forces

for n = 1:nnbc
    nb = inbc(n);
    g_f(nb) = g_f(nb)+vnbc(n);
end

% Apply known displacements (1-0) method

for nj = 1:nebc
    j = iebc(nj);
    for k=1:neq
        if k~=j
            g_f(k) = g_f(k)-g_k(k,j)*vebc(nj);
            g_k(k,j) = 0;
            g_k(j,k) = 0;
        else
            g_k(j,j) = 1;
            g_f(j) = vebc(nj);
        end
    end
end
sol = g_k\g_f;

% to check solutions

TT = table(x,sol);
