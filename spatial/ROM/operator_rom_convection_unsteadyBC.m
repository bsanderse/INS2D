function [conv_bc,conv_linear,conv_quad] = operator_rom_convection_unsteadyBC(P,options)

Cux = options.discretization.Cux;
Cuy = options.discretization.Cuy;
Cvx = options.discretization.Cvx;
Cvy = options.discretization.Cvy;

Au_ux = options.discretization.Au_ux;
Au_uy = options.discretization.Au_uy;
Av_vx = options.discretization.Av_vx;
Av_vy = options.discretization.Av_vy;

yAu_ux = options.discretization.yAu_ux;
yAu_uy = options.discretization.yAu_uy;
yAv_vx = options.discretization.yAv_vx;
yAv_vy = options.discretization.yAv_vy;

Iu_ux = options.discretization.Iu_ux;
Iv_uy = options.discretization.Iv_uy;
Iu_vx = options.discretization.Iu_vx;
Iv_vy = options.discretization.Iv_vy;

yIu_ux = options.discretization.yIu_ux;
yIv_uy = options.discretization.yIv_uy;
yIu_vx = options.discretization.yIu_vx;
yIv_vy = options.discretization.yIv_vy;

A = blkdiag([Au_ux; Au_uy],[Av_vx; Av_vy]); %not efficient but intuitive
I = [blkdiag(Iu_ux,Iv_uy);blkdiag(Iu_vx,Iv_vy)]; %not efficient but intuitive
K = blkdiag([Cux Cuy],[Cvx Cvy]); %not efficient but intuitive
y_A = [yAu_ux;yAu_uy;yAv_vx;yAv_vy];
y_I = [yIu_ux;yIv_uy;yIu_vx;yIv_vy];



