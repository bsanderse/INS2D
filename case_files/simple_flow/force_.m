function [Fx,Fy,dFx,dFy] = force_(~,t,options,~)

% Nux_in = options.grid.Nux_in;
% Nuy_in = options.grid.Nuy_in;
% Nvx_in = options.grid.Nvx_in;
% Nvy_in = options.grid.Nvy_in;
% 
% Fx = zeros(Nux_in,Nuy_in);
% Fy = zeros(Nvx_in,Nvy_in);
% dFx = zeros(Nux_in,Nuy_in);
% dFy = zeros(Nvx_in,Nvy_in);

Nu = options.grid.Nu;
Nv = options.grid.Nv;

Fx = zeros(Nu,1);
Fy = zeros(Nv,1);
dFx = zeros(Nu,1);
dFy = zeros(Nv,1);