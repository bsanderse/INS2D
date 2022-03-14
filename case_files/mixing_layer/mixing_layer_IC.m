function [u,v,p,options] = mixing_layer_IC(t,options)
% initial velocity field actuator

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

% extend inflow 
Delta_U = 1;
u = 1 + (Delta_U/2)*tanh(2*options.grid.yu);
% u  = ones(Nux_in,Nuy_in);  
v  = zeros(Nvx_in,Nvy_in);

% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
    
end