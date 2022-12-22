function [u,v,p,T,options] = RBC_EC_IC(t,options)
% initial velocity field RBC: load state state results

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

% load velocity field and temperature
results = load('results/RBC_Ra1e4_Ge_1e-1_32x32/steady_state.mat');
u  = results.V(options.grid.indu); %zeros(Nux_in,Nuy_in);
v  = results.V(options.grid.indv); %zeros(Nvx_in,Nvy_in);
T  = results.T;
% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p  = zeros(Npx,Npy);
 
% temperature: linear profile from 1 (bottom) to 0 (top)
% T  = (Ly - ypp)/Ly;
% inverted temperature: linear profile from 1 (bottom) to 0 (top)
% T  = 0.5 - (Ly - ypp)/Ly ;
% homogeneous 0.5
% T  = 0.5*ones(Npx,Npy);
% homogeneous 0
% T  = zeros(Npx,Npy); 
% high order polynomial that has zero derivatives
%T = 10*(ypp.^2) .* ((Ly - ypp)/Ly).^2;
% T = 2*ypp.^3-3*ypp.^2 + 1;
% polynomial that is "periodic"
% T = ypp.^2 .* (1 - ypp).^2;
% 
% T = sin(2*pi*ypp); %
end