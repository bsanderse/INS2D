function [maxres,Fres,dF] = F(V,C,p,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
% V: velocity field
% C: 'convection' field: e.g. d(c_x u)/dx + d(c_y u)/dy; usually c_x = u,
% c_y=v
% p: pressure

if (nargin<6)
    getJacobian = 0;
end

Nu = options.grid.Nu;
Nv = options.grid.Nv;

% unsteady BC
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end

% pressure
Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
y_px = options.discretization.y_px;
y_py = options.discretization.y_py;
Gpx  = Gx*p + y_px;
Gpy  = Gy*p + y_py;

% convection:
[convu, convv, dconvu, dconvv] = convection(V,C,t,options,getJacobian);

% diffusion
[d2u, d2v, dDiffu, dDiffv] = diffusion(V,t,options,getJacobian);

% body force
[Fx, Fy] = force(t,options);

% residual in Finite Volume form, including the pressure contribution
Fu   = - convu + d2u - Gpx + Fx;
Fv   = - convv + d2v - Gpy + Fy;

Fres = [Fu;Fv];

% norm of residual
maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
        
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
    
    dFu  = - dconvu + dDiffu;
    dFv  = - dconvv + dDiffv;
    
    dF   = [dFu; dFv];
    
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end