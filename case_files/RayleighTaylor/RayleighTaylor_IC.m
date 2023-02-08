function [u,v,p,T,options] = RayleighTaylor_IC(t,options)
% initial velocity field RBC

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
% xpp = options.grid.xpp;
ypp = options.grid.ypp;
% Lx = options.grid.x2 - options.grid.x1;
Ly = options.grid.y2 - options.grid.y1;

% constant velocity field
u  = zeros(Nux_in,Nuy_in);
v  = zeros(Nvx_in,Nvy_in);

% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
 
% inverted profile for temperature
T  = ypp<(Ly/2);

end