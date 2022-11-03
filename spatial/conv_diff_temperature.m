function [conv_diff, Jac_conv_diff_T, Jac_conv_diff_V] = conv_diff_temperature(T,V,t,options,getJacobian)
% evaluate convection and diffusion terms for temperature equation and optionally Jacobian
% the Jacobian is with respect to both T and V
% note that all terms are in integrated (finite volume) form

% visc = options.case.visc;

indu = options.grid.indu;
indv = options.grid.indv;

uh = V(indu);
vh = V(indv);

% Nu = options.grid.Nu;
% Nv = options.grid.Nv;
NV = options.grid.NV;
NT = options.grid.NT;

Jac_conv_diff_T = spalloc(NT,NT,0);
Jac_conv_diff_V = spalloc(NT,NV,0);

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

u_Tx = Iu_Tx*uh + yIu_Tx;
T_Tx = AT_Tx*T + yAT_Tx;
v_Ty = Iv_Ty*vh + yIv_Ty;
T_Ty = AT_Ty*T + yAT_Ty;
conv_diff = -CTx*(u_Tx.*T_Tx) - CTy*(v_Ty.*T_Ty) + ...
            (DiffT*T + yDiffT);

        
if (getJacobian == 1)
    N1 = length(u_Tx); %options.grid.N1;
    N2 = length(v_Ty); %options.grid.N2;
    N3 = length(T_Tx); %options.grid.N3;
    N4 = length(T_Ty); %options.grid.N4;
    
    % convective terms, d/dx (u*T)
    C1         = CTx*spdiags(u_Tx,0,N1,N1);
    Conv_Tx    = C1*AT_Tx;

    % convective terms, d/dy (v*T)
    C2         = CTy*spdiags(v_Ty,0,N2,N2);
    Conv_Ty    = C2*AT_Ty;

    % diffusive terms: Jacobian is simply DiffT
    % Jacobian with respect to T is then:
    Jac_conv_diff_T = -(Conv_Tx + Conv_Ty) + DiffT;

    % Jacobian with respect to V:
    % convective terms, d/dx (u*T)
    C1         = CTx*spdiags(T_Tx,0,N3,N3);
    Conv_u     = C1*Iu_Tx;

    % convective terms, d/dy (v*T)
    C2         = CTy*spdiags(T_Ty,0,N4,N4);
    Conv_v     = C2*Iv_Ty;    

    Jac_conv_diff_V = -[Conv_u Conv_v];

else



end

    

end

