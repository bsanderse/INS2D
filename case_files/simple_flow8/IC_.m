function [u,v,p,options] = IC_(t,options)
% initial velocity field actuator

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

% actuator, extend inflow 
% u  = ones(Nux_in,Nuy_in);  %wrong for manipulated BC
% v  = zeros(Nvx_in,Nvy_in);

% options = set_bc_vectors(0,options);
% f       = options.discretization.yM;
% dp      = pressure_poisson(f,0,options);
% V0      = - options.grid.Om_inv.*(options.discretization.G*dp);
% u = V0(1:Nux_in*Nuy_in);
% v = V0(Nux_in*Nuy_in+1:end);

u = ones(Nux_in*Nuy_in,1);
% u = zeros(Nux_in*Nuy_in,1);
v = zeros(Nvx_in*Nvy_in,1);

u = reshape(u,Nux_in,Nuy_in);
v = reshape(v,Nvx_in,Nvy_in);

% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
    
end