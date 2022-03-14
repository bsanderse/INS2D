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

Om = options.grid.Om;
Om_u = Om(1:Nu);
Om_v = Om(Nu+1:Nv);

% Fx = zeros(Nu,1);
Fy = zeros(Nv,1);

Fx = Om_u.*ones(Nu,1);
% Fy = ones(Nv,1);

dFx = zeros(Nu,1);
dFy = zeros(Nv,1);