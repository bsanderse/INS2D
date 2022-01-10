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
K_h = options.discretization.K_h;
I_h = options.discretization.I_h;
A_h = options.discretization.A_h;
y_I = options.discretization.y_I;
id_normal = options.grid.id_normal;
id_tangential = options.grid.id_tangential;
id_n_t = id_normal+id_tangential;
V_n_t = id_n_t.*V;
gO = options.BC.gO;
% gO = @(u) 0; % botch

% % conv_Diag1 = diag(K_h*diag(I_h*V+y_I)*A_h);
% conv_Diag11 = diag(K_h*((I_h*V+y_I).*A_h));
NF = length(y_I);
conv_Diag1 = diag(K_h*spdiags(I_h*V+y_I,0,NF,NF)*A_h);
% conv_Diag12 = diag(K_h*(spdiags(I_h*V+y_I,0,NF,NF)*A_h));
% norm(conv_Diag12-conv_Diag1)
% norm(conv_Diag11-conv_Diag1)
% % conv_Diag1 = dot(K_h',(I_h*V+y_I).*A_h)';

%% performance testing
% profile on
% (I_h*V+y_I).*A_h;
% spdiags(I_h*V+y_I,0,NF,NF)*A_h;
% 
% (I_h*V+y_I).*A_h(:,1);
% spdiags(I_h*V+y_I,0,NF,NF)*A_h(:,1);
% 
% (y_I).*A_h(:,1);
% spdiags(y_I,0,NF,NF)*A_h(:,1);
% 
% profile off
% profile viewer
%%

y_O1 = conv_Diag1.*V_n_t; % too (storage+time) expensive implementation
% y_O2 = diag(K_h*((I_h*V+y_I).*A_h)).*V_n_t; % too (time) expensive implementation
% y_O1 = y_O2;
% y_O3 = dot(K_h',(I_h*V+y_I).*A_h)'.*V_n_t; % too (time) expensive implementation
Conv_diag = options.grid.C;
y_O_diag = (Conv_diag*V).*V;

y_O  = y_O_diag - V_n_t.*gO(V);
Fres = Fres + y_O;
%% tests
Conv = [convu; convv];
if (V'*(-Conv + y_O_diag )>1e-14)
    warning('convection skew-symmetry not fixed by mvp-obc') %only works in absence of Dirichlet BC
    V'*(-Conv + y_O_diag)
end
% V'*(-Conv + y_O_diag )
% norm(y_O1-y_O_diag)  %error of 1e-5 for 3x3 grid, error of 1e-7 for 20x20
% norm(conv_Diag1.*V-y_O_diag)
%%
options.solversettings.Newton_factor = 3;
[~,~,~,~, diag_] = convection(V,V,t,options,1);
% profile off
% profile viewer

hx = options.grid.hx;
NF = length(y_I);


% norm of residual
maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    % we return only the Jacobian with respect to V (not p)
           
    dFu  = - dconvu + dDiffu + dFx;
    dFv  = - dconvv + dDiffv + dFy;
    
    dF   = [dFu; dFv];
    
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end