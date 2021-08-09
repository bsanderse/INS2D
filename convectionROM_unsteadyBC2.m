function [Conv, Jac] = convectionROM_unsteadyBC2(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians
Jac = -666;
if getJacobian
    error('Jacobian not implemented in offline decomposition')
end

% order4     = options.discretization.order4;
% regularize = options.case.regularize;
% 
% M   = options.rom.M;
% B = options.rom.B;

% V = getFOM_velocity(R,t,options);
% C = V;

%% Henrik's pfusch
% Cux = options.discretization.Cux;
% Cuy = options.discretization.Cuy;
% Cvx = options.discretization.Cvx;
% Cvy = options.discretization.Cvy;
% 
% Au_ux = options.discretization.Au_ux;
% Au_uy = options.discretization.Au_uy;
% Av_vx = options.discretization.Av_vx;
% Av_vy = options.discretization.Av_vy;
% 
% yAu_ux = options.discretization.yAu_ux;
% yAu_uy = options.discretization.yAu_uy;
% yAv_vx = options.discretization.yAv_vx;
% yAv_vy = options.discretization.yAv_vy;
% 
% Iu_ux = options.discretization.Iu_ux;
% Iv_uy = options.discretization.Iv_uy;
% Iu_vx = options.discretization.Iu_vx;
% Iv_vy = options.discretization.Iv_vy;
% 
% yIu_ux = options.discretization.yIu_ux;
% yIv_uy = options.discretization.yIv_uy;
% yIu_vx = options.discretization.yIu_vx;
% yIv_vy = options.discretization.yIv_vy;

% A = blkdiag([Au_ux; Au_uy],[Av_vx; Av_vy]); %not efficient but intuitive
% I = [blkdiag(Iu_ux,Iv_uy);blkdiag(Iu_vx,Iv_vy)]; %not efficient but intuitive
% K = blkdiag([Cux Cuy],[Cvx Cvy]); %not efficient but intuitive
% y_A = [yAu_ux;yAu_uy;yAv_vx;yAv_vy];
% y_I = [yIu_ux;yIv_uy;yIu_vx;yIv_vy];

% Conv = B'*K*(((I*C)+y_I).*((A*V)+y_A));

% B = options.rom.B;
% Conv_hom =  B'*K*((I*B*R).*(A*B*R));
% Conv_hom2 = options.rom.C_hom*kron(R,R);
% norm(Conv_hom-Conv_hom2)

R_inhom = get_a_inhom(t,options);
R_bc    = get_a_bc(t,options);

Conv = options.rom.C_hom*kron(R,R) ...
     + options.rom.C_hom_inhom*kron(R,R_inhom) ...
     + options.rom.C_hom_bc*kron(R,R_bc) ...
     + options.rom.C_inhom*kron(R_inhom,R_inhom) ...
     + options.rom.C_inhom_bc*kron(R_inhom,R_bc) ...
     + options.rom.C_bc*kron(R_bc,R_bc);