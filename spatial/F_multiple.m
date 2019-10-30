function [maxres,Fres,dF] = F_multiple(V,C,p,t,options,getJacobian)
% call F.m for multiple (V,p) pairs, as required for example in implicit RK
% methods

if (nargin<5)
    getJacobian = 0;
end

s_RK = options.time.s_RK;
% c_RK = options.time.c_RK;

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu + Nv;
Np = options.grid.Np;

maxres = zeros(s_RK,1);
Fres   = zeros(s_RK*NV,1);
dF     = [];

for i=1:s_RK
       
    % take stage i out of total vector:
    indxV = (1:NV) + NV*(i-1);
    indxp = (1:Np) + Np*(i-1);
    Vi    = V(indxV);
    Ci    = C(indxV);
    pi    = p(indxp);
    ti    = t(i);
    
    % compute residual and Jacobian for this stage
    [maxresi,Fresi,dFi] = F(Vi,Ci,pi,ti,options,getJacobian);

    maxres(i)   = maxresi;
    Fres(indxV) = Fresi;
    % add Jacobian in block diagonal form; this allows multiplication with A_RK
    % later on
    dF          = blkdiag(dF,dFi);
end