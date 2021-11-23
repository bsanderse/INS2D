function [maxres,Fres,dFres] = F_ROM_notvelocityonly(R,p,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field

% warning pressure contribution for outflow BC not included
if (nargin<5)
    getJacobian = 0;
end

if (getJacobian == 1)
    error('bc recon = 2 and get Jacobian = 1 not implemented')
end
V = getFOM_velocity(R,t,options);
p = pressure_additional_solve(V,p,t,options);
nopressure = 0;

[maxres,Fres,dFres] = myF(V,V,p,t,options,getJacobian,nopressure);
B = options.rom.B;
Fres = B'*Fres;

maxres  = max(abs(Fres));


end