function [conv_diff, Jac_conv_diff] = conv_diff_temperature(T,V,t,options,getJacobian)
% evaluate convection and diffusion terms for temperature equation and optionally Jacobian
% note that all terms are in integrated (finite volume) form

% visc = options.case.visc;

indu = options.grid.indu;
indv = options.grid.indv;

uh = V(indu);
vh = V(indv);

NT = options.grid.NT;

Jac_conv_diff = spalloc(NT,NT,0);

CTx = options.discretization.CTx;
CTy = options.discretization.CTy;

AT_Tx = options.discretization.AT_Tx;
AT_Ty = options.discretization.AT_Ty;
yAT_Tx = options.discretization.yAT_Tx;
yAT_Ty = options.discretization.yAT_Ty;


Iu_Tx = options.discretization.Iu_Tx;
Iv_Ty = options.discretization.Iv_Ty;
yIu_Tx = options.discretization.yIu_Tx;
yIv_Ty = options.discretization.yIv_Ty;
 
DiffT = options.discretization.DiffT;
yDiffT = options.discretization.yDiffT;

        
conv_diff = -CTx*((Iu_Tx*uh + yIu_Tx).*(AT_Tx*T + yAT_Tx)) + ...
            -CTy*((Iv_Ty*vh + yIv_Ty).*(AT_Ty*T + yAT_Ty)) + ...
            (DiffT*T + yDiffT);

        
if (getJacobian == 1)
    
    
end

    

end

