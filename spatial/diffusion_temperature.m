function [diff, Jac_diff_T, Jac_diff_V] = diffusion_temperature(T,t,options,getJacobian)
% evaluate convection and diffusion terms for temperature equation and optionally Jacobian
% the Jacobian is with respect to both T and V
% note that all terms are in integrated (finite volume) form

% visc = options.case.visc;

% indu = options.grid.indu;
% indv = options.grid.indv;

% Nu = options.grid.Nu;
% Nv = options.grid.Nv;
NV = options.grid.NV;
NT = options.grid.NT;

Jac_diff_T = spalloc(NT,NT,0);
Jac_diff_V = spalloc(NT,NV,0);


DiffT  = options.discretization.DiffT;
yDiffT = options.discretization.yDiffT;

% note that DiffT and yDiffT include already the scaling with alfa4
diff   = (DiffT*T + yDiffT);

        
if (getJacobian == 1)
    % diffusive terms: Jacobian is simply DiffT
    % Jacobian with respect to T is then:
    Jac_diff_T = DiffT;

end

    

end

