function [maxres,Fu,Fv] = F(uh,vh,p,t,options)
% calculate rhs of momentum equations 

Omu_inv = options.grid.Omu_inv;
Omv_inv = options.grid.Omv_inv;

Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
y_px = options.discretization.y_px;
y_py = options.discretization.y_py;

% convection:
[convu, convv] = convection(uh,vh,t,options);

% diffusion
[d2u, d2v] = diffusion(uh,vh,t,options);

% body force
[Fx, Fy] = force(t,options);

% residual
Fu   = Omu_inv.*(- convu + d2u - Gx*p - y_px + Fx);
Fv   = Omv_inv.*(- convv + d2v - Gy*p - y_py + Fy);

maxres  = max(abs([Fu;Fv]));

end