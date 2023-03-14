function [u,v,p,T,options] = RBC_steady_IC(t,options)
% initial velocity field RBC

Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
% xpp = options.grid.xpp;
ypp = options.grid.ypp;
xu = options.grid.xu;
yu = options.grid.yu;
xv = options.grid.xv;
yv = options.grid.yv;
% Lx = options.grid.x2 - options.grid.x1;
Ly = options.grid.y2 - options.grid.y1;

% constant velocity field
% u  = zeros(Nux_in,Nuy_in);
% v  = zeros(Nvx_in,Nvy_in);

load_IC = 1;

if (load_IC)
    
    % load IC at different Gebhart number
    filename = ['results/RBC/runs_Mar2023/T_Ra1e5_Ge1e-1_N' num2str(Npx) '.mat'];
    disp(['loading ' filename]);
    data_IC = load(filename);
    u = data_IC.V(options.grid.indu);
    v = data_IC.V(options.grid.indv);
    p = data_IC.p;
    T = data_IC.T;
    
    
else
    
    % constant velocity field
    % u  = zeros(Nux_in,Nuy_in);
    % v  = zeros(Nvx_in,Nvy_in);


    % regularized lid-driven cavity field, see T.M. Shih, C.H. Tan, and B.C. Hwang. Effects of grid staggering on numerical schemes.
    % International Journal for Numerical Methods in Fluids, 9:192–212, 1989.
    u  = -8*(xu.^4 - 2*xu.^3 + xu.^2).*(4*yu.^3 - 2*yu);
    v  = +8*(4*xv.^3 - 6*xv.^2 + 2*xv).*(yv.^4 - yv.^2);

    % pressure: should in principle NOT be prescribed. will be calculated if
    % p_initial=1
    p  = zeros(Npx,Npy);

    % T = load
    % temperature: linear profile from 1 (bottom) to 0 (top)
%     T  = (Ly - ypp)/Ly;
    % inverted temperature: linear profile from 1 (bottom) to 0 (top)
%     T  = 0.5 - (Ly - ypp)/Ly ;
    % homogeneous 0.5
    % T  = 0.5*ones(Npx,Npy);
    % homogeneous 0
    % T  = zeros(Npx,Npy); 
    % high order polynomial that has zero derivatives
    %T = 10*(ypp.^2) .* ((Ly - ypp)/Ly).^2;
    % T = 2*ypp.^3-3*ypp.^2 + 1;
    % T = ypp.^5;
    % polynomial that is "periodic"
    % T = ypp.^2 .* (1 - ypp).^2;
    % 
    % T = sin(2*pi*ypp); %
    % random between 0 and 1
    T  = rand(Npx,Npy);
    
end
end