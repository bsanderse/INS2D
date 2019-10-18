function [u,v,p,options] = actuator_ROM_IC(t,options)
% initial velocity field actuator

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

yu = options.grid.yu;
% actuator, extend inflow 
% u  = ones(Nux_in,Nuy_in);  
v  = zeros(Nvx_in,Nvy_in);

u  = 0.75-(3/32)*(yu-2).*(yu+2);  


% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
    
end