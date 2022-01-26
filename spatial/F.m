function [maxres,Fres,dF] = F(V,C,p,t,options,getJacobian,nopressure)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
% V: velocity field
% C: 'convection' field: e.g. d(c_x u)/dx + d(c_y u)/dy; usually c_x = u,
% c_y=v, so C=V
% p: pressure

% getJacobian = 1: return dFdV
% nopressure = 1: exclude pressure gradient; in this case input argument p is not used

if (nargin<7)
    nopressure = 0;    
end
if (nargin<6)
    getJacobian = 0;
end

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = options.grid.NV;

% unsteady BC
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end

if (nopressure == 0)
    % pressure
    Gx   = options.discretization.Gx;
    Gy   = options.discretization.Gy;
    y_px = options.discretization.y_px;
    y_py = options.discretization.y_py;
    Gpx  = Gx*p + y_px;
    Gpy  = Gy*p + y_py;
end

% convection:
[convu, convv, dconvu, dconvv] = convection(V,C,t,options,getJacobian);

% diffusion
% [d2u, d2v, dDiffu, dDiffv] = diffusion(V,t,options,getJacobian);
[d2u, d2v, dDiffu, dDiffv] = mydiffusion(V,t,options,getJacobian);

% body force
if (options.force.isforce == 1)
    if (options.force.force_unsteady == 1)
        [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian);
    else
        Fx = options.force.Fx;
        Fy = options.force.Fy;
        dFx = spalloc(Nu,NV,0);
        dFy = spalloc(Nv,NV,0);            
    end
else
    Fx  = zeros(Nu,1);
    Fy  = zeros(Nv,1);
    dFx = spalloc(Nu,NV,0);
    dFy = spalloc(Nv,NV,0);      
end

% residual in Finite Volume form, including the pressure contribution
Fu   = - convu + d2u + Fx;
Fv   = - convv + d2v + Fy;

% nopressure=0 is the most common situation, in which we return the entire 
% right-hand side vector
if (nopressure == 0) 
    Fu = Fu - Gpx;
    Fv = Fv - Gpy;
end

Fres = [Fu;Fv];

%% mvp-obc 
%%% problem: fixes skew-symmetry for all diagonal entries even without obc!
%%% -> solved when we use options.grid.C
% K_h = options.discretization.K_h;
% I_h = options.discretization.I_h;
% A_h = options.discretization.A_h;
% y_I = options.discretization.y_I;
id_normal = options.grid.id_normal;
id_tangential = options.grid.id_tangential;
id_n_t = id_normal+id_tangential;
V_n_t = id_n_t.*V;
gO = options.BC.gO;
dgO = options.BC.dgO;
% gO = @(u) 0; % botch

% NF = length(y_I);
% % conv_matrix = K_h*spdiags(I_h*V+y_I,0,NF,NF)*A_h;
% conv_Diag1 = diag(K_h*spdiags(I_h*V+y_I,0,NF,NF)*A_h);

Conv_diag = options.grid.C;
y_O_diag = Conv_diag*V;

y_O  = y_O_diag.*V - V_n_t.*gO(V);
% y_O  = conv_Diag1.*V- V_n_t.*gO(V);
Fres = Fres + y_O;

%% tests
% Conv = [convu; convv];
% if (V'*(-Conv + y_O_diag.*V )>1e-14)
%     warning('convection skew-symmetry not fixed by mvp-obc') %only works in absence of Dirichlet BC
%     V'*(-Conv + y_O_diag.*V)
% end

% if (V'*(-Conv + conv_Diag1.*V )>1e-14)
%     warning('convection skew-symmetry not fixed by mvp-obc') %only works in absence of Dirichlet BC
%     V'*(-Conv + conv_Diag1.*V)
% end

% V'*(-Conv + y_O_diag.*V )
% V'*(-Conv + conv_Diag1.*V )
% norm(conv_Diag1-y_O_diag)
% 
% y_M = options.discretization.yM;
% M_h = options.discretization.M;
% norm(M_h*V+y_M)
% 
% y_A = options.discretization.y_A;
% y_conv = K_h*spdiags(I_h*V+y_I,0,NF,NF)*y_A;
% Conv2 = conv_matrix*V + y_conv;
% norm(Conv-Conv2)
% 
% non_diag = conv_matrix-diag(diag(conv_matrix));
% V'*non_diag*V
% norm(full(non_diag+non_diag'))
% 
% hx = options.grid.hx;
% full([conv_Diag1-y_O_diag conv_Diag1 y_O_diag V*hx(1)/2
% ])

% elf = 11;
% pro_K11 = K_h(elf,:);
% vec = find(K11);
% K11 = pro_K11(vec);
% I11 = I_h(vec,:);
% A11 = A_h(vec,elf);
% V11 = V(vec);
% K11*((I11*V).*A11)
%%

% norm of residual
maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    % we return only the Jacobian with respect to V (not p)
               
    dFu  = - dconvu + dDiffu + dFx;
    dFv  = - dconvv + dDiffv + dFy;
    
    dF   = [dFu; dFv];

    Jac_y_O = diag(Conv_diag*V) + diag(V)*Conv_diag ...
            - diag(id_n_t)*diag(gO(V)) - diag(id_n_t)*diag(dgO(V))*diag(V);
    dF   = dF + Jac_y_O;
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end