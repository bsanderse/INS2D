function [maxres,Fres,dF] = F(V,p,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
if (nargin<5)
    getJacobian = 0;
end

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu+Nv;

uh = V(1:Nu);
vh = V(Nu+1:Nu+Nv);

Om_inv = options.grid.Om_inv;

Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
y_px = options.discretization.y_px;
y_py = options.discretization.y_py;

% convection:
[convu, convv, dconvu, dconvv] = convection(uh,vh,t,options,getJacobian);

% diffusion
[d2u, d2v, dDiffu, dDiffv] = diffusion(uh,vh,t,options,getJacobian);

% body force
[Fx, Fy] = force(t,options);

% residual
Fu   = - convu + d2u - Gx*p - y_px + Fx;
Fv   = - convv + d2v - Gy*p - y_py + Fy;

if (options.case.steady==0) % unsteady case, solve for velocities
    Fres = Om_inv.*[Fu;Fv];
else
    Fres = [Fu;Fv];
end

maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
        
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
    
    dFu  = - dconvu + dDiffu;
    dFv  = - dconvv + dDiffv;
    
    dF   = [dFu; dFv];
    
    if (options.case.steady==0) % unsteady case, solve for velocities
        dF = spdiags(Om_inv,0,NV,NV)*dF;
    end
    
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end