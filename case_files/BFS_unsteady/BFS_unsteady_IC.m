function [u,v,p,options] = BFS_unsteady_IC(t,options)
% initial velocity field LDC

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

yu = options.grid.yu;

% BFS, extend inflow 
u  = (yu>=0).*(24*yu.*(1/2-yu));  
v  = zeros(Nvx_in,Nvy_in);

% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
    
end