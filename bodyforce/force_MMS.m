%% MMS of Eca for k-epsilon

% discrete differentiation of exact solution:
% u   = erf(sigma*yu./xu);
% v   = 1/(sigma*sqrt(pi)) * (1 - exp(-(sigma*yv./xv).^2));
% pres = 0.5*log(2*xpp-xpp.^2+0.25).*log(4*ypp.^3-3*ypp.^2+1.25);
% % pres= 0.5*log(0.5*xpp-0.75*xpp.^2+(1/3)*xpp.^3+0.25).*...
% %           log(4*ypp.^3-3*ypp.^2+1.25);
% kt  = kmax*(sigma_nu*ypp./xpp).^2 .* exp(1-(sigma_nu*ypp./xpp).^2);
% e   = 0.36*kmax^2/numax*exp(-(sigma_nu*ypp./xpp).^2);
% 
% uh  = u(:);  
% vh  = v(:);
% V   = [uh;vh];
% p   = pres(:);
% kth = kt(:);
% eh  = e(:);
% 
% % convection of momentum
% fuconv =   Cux*((Iu_ux*uh+yIu_ux).*(Au_ux*uh+yAu_ux)) ...
%          + Cuy*((Iv_uy*vh+yIv_uy).*(Au_uy*uh+yAu_uy));
% fvconv =   Cvx*((Iu_vx*uh+yIu_vx).*(Av_vx*vh+yAv_vx)) ...
%          + Cvy*((Iv_vy*vh+yIv_vy).*(Av_vy*vh+yAv_vy));
% 
% % diffusion of momentum
% nu_ux      = spdiags( nu + Cmu* (Ak_ux*kth+yAk_ux).^2 ./ (Ae_ux*eh+yAe_ux + eps),0,N1,N1);
% nu_uy      = spdiags( nu + Cmu* (Ak_uy*kth+yAk_uy).^2 ./ (Ae_uy*eh+yAe_uy + eps),0,N2,N2);
% nu_vx      = spdiags( nu + Cmu* (Ak_vx*kth+yAk_vx).^2 ./ (Ae_vx*eh+yAe_vx + eps),0,N3,N3);
% nu_vy      = spdiags( nu + Cmu* (Ak_vy*kth+yAk_vy).^2 ./ (Ae_vy*eh+yAe_vy + eps),0,N4,N4);
%             
% fudiff     = Dux*( nu_ux* 2 * (Su_ux*uh + ySu_ux)) + ...
%              Duy*( nu_uy* (Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy));
% fvdiff     = Dvx*( nu_vx* (Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx)) + ...
%              Dvy*( nu_vy* 2 * (Sv_vy*vh + ySv_vy));
% 
% % pressure terms
% fupres = Gx*p + y_px;
% fvpres = Gy*p + y_py;
% 
% Fx     = fuconv + fupres - fudiff;
% Fy     = fvconv + fvpres - fvdiff;
% 
% 
% % source terms for k
% % convection of k
% fkconv = Ckx*( (Iu_kx*uh+yIu_kx).*(Ak_kx*kth+yAk_kx) ) + ...
%          Cky*( (Iv_ky*vh+yIv_ky).*(Ak_ky*kth+yAk_ky) );
% 
% % diffusion of k
% nu_kx  = Cmu * (Ak_kx*kth + yAk_kx).^2 ./(Ae_ex*eh+yAe_ex + eps); 
% nu_ky  = Cmu * (Ak_ky*kth + yAk_ky).^2 ./(Ae_ey*eh+yAe_ey + eps);
% fkdiff = Dkx*( ( nu + nu_kx/sigmak ) .* (Skx*kth+ySkx) ) + ...
%          Dky*( ( nu + nu_ky/sigmak ) .* (Sky*kth+ySky) ); 
% 
% % production of k
% fkprod = Omp.*Cmu.*(kth.^2 ./eh) .* ( ...
%            2*(Cux_k*uh+yCux_k).^2 + 2*(Cvy_k*vh+yCvy_k).^2 + ...
%              (Cuy_k*(Auy_k*uh+yAuy_k)+yCuy_k + Cvx_k*(Avx_k*vh+yAvx_k)+yCvx_k).^2 );
% 
% % dissipation of k
% fkdiss = Omp.*eh;
% 
% Fk     = fkconv - fkdiff - fkprod + fkdiss;
% 
% % source terms for e
% % convection of e
% feconv = Ckx*( (Iu_kx*uh+yIu_kx).*(Ae_ex*eh+yAe_ex) ) + ...
%          Cky*( (Iv_ky*vh+yIv_ky).*(Ae_ey*eh+yAe_ey) );
% 
% % diffusion of e
% fediff = Dkx*( ( nu + nu_kx/sigmae ) .* (Sex*eh+ySex) ) + ...
%          Dky*( ( nu + nu_ky/sigmae ) .* (Sey*eh+ySey) ); 
%      
% % production of e
% feprod = Ce1*eh./kth.*fkprod;
% 
% % dissipation of e
% fediss = Ce2*eh./kth.*fkdiss;
% 
% Fe     = feconv - fediff - feprod + fediss;