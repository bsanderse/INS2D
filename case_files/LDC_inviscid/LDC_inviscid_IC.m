function [u,v,p,options] = LDC_inviscid_IC(t,options)
    % initial velocity field LDC

    Npx = options.grid.Npx;
    Npy = options.grid.Npy;
    Nux_in = options.grid.Nux_in;
    Nuy_in = options.grid.Nuy_in;
    Nvx_in = options.grid.Nvx_in;
    Nvy_in = options.grid.Nvy_in;

    % random velocity field
    u  = -1 + 2*rand(Nux_in,Nuy_in);
    v  = -1 + 2*rand(Nvx_in,Nvy_in);
%     u  = zeros(Nux_in,Nuy_in);
%     v  = zeros(Nvx_in,Nvy_in);

    % pressure: should in principle NOT be prescribed. will be calculated if
    % p_initial=1
    p  = zeros(Npx,Npy);

end