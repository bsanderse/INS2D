function [u,v] = LDC_IC(t,options)
% initial velocity field LDC

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

% constant:    
u  = ones(Nux_in,Nuy_in);
v  = zeros(Nvx_in,Nvy_in);

    
end