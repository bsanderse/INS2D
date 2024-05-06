function [u,v,p,options] = IC_(t,options)
% initial velocity field actuator

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

xu = options.grid.yu;

A = options.case.IC_params(1);
f   = options.case.IC_params(2);
phi   = options.case.IC_params(3);

x0_ = @(omega,A,f,phi) A*sin(2*pi*f*omega + phi);
% As = [.8, .9, 1., 1.1, 1.2];
% fs = [1 2 3];
% phis = [-.25, -.125, 0 , .125 .25];

u = x0_(xu,A,f,phi);
v = zeros(Nvx_in,Nvy_in);


% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
    
end