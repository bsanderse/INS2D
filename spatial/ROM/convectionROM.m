function [Conv, Jac] = convectionROM(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

order4     = options.discretization.order4;
regularize = options.case.regularize;

M   = options.rom.M;

V = getFOM_velocity(R,t,options);
C = V;

%% Henrik's pfusch
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

Conv = K*(((I*C)+y_I).*((A*V)+y_A));


% if (order4==0)
%     
%     %% no regularization
%     if (regularize == 0)
%         
%         Conv_quad   = options.rom.Conv_quad;
%         Conv_linear = options.rom.Conv_linear;
%         yConv       = options.rom.yConv;
%         
%         Conv        = Conv_quad*kron(R,R) + Conv_linear*R + yConv;
%          
%         if (getJacobian == 1)          
%            % d/dR (kron(R,R)) = kron(E,R) + kron(R,E), where E = speye(M)
%            E   = speye(M);
%            Jac = Conv_quad*(kron(E,R) + kron(R,E)) + Conv_linear;
%         else
%            Jac = 0;
% %            Jac = spalloc(M,M,0);
%         end
%         
%     else
%         error('not implemented');
%     end
%  
% 
% else
%     
%     error('not implemented');
%     
% end