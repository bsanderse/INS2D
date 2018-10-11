function [maxres,Fres,dF] = F_multiple(V,p,t,options,getJacobian)
% call F.m for multiple (V,p) pairs, as required for example in implicit RK
% methods

if (nargin<5)
    getJacobian = 0;
end

s_RK = options.time.s_RK;
c_RK = options.time.c_RK;

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu + Nv;
Np = options.grid.Np;

maxres = [];
Fres   = [];
dF     = [];

for i=1:s_RK
    
%     indxV = (1:NV) + (NV+Np)*(i-1);
%     indxp = (NV+1:NV+Np) + (NV+Np)*(i-1);
    indxV = (1:NV) + NV*(i-1);
    indxp = (1:Np) + Np*(i-1);
    Vi = V(indxV);
    pi = p(indxp);
    ti = c_RK(i);
    [maxresi,Fresi,dFi] = F(Vi,pi,ti,options,getJacobian);

    maxres = [maxres; maxresi];
    Fres   = [Fres; Fresi];
    dF     = blkdiag(dF,dFi);
end