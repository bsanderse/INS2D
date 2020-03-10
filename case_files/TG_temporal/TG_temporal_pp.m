%% post-processing Taylor-Green
xu = options.grid.xu;
yu = options.grid.yu;
xv = options.grid.xv;
yv = options.grid.yv;
xpp = options.grid.xpp;
ypp = options.grid.ypp;

Nu = options.grid.Nu;
Nv = options.grid.Nv;
Np = options.grid.Np;

uh = V(options.grid.indu);
vh = V(options.grid.indv);

% temporal = 0;
% % keyboard;
% if (temporal==1)
%     load('results/Taylor_Green/TG_Re100_N20_dt1e-3_Dir.mat');
% %     load('results/Taylor_Green/TG_Re100_N20_dt1e-3_per.mat');
%     u_error = uh-ufine;
%     v_error = vh-vfine;
%     
%     p       = p-mean(p(:));
%     p_error = p-pfine;
%     
%     u_error_i(jj,j) = max(abs(u_error))
%     v_error_i(jj,j) = max(abs(v_error))
%     u_error_2(jj,j) = sqrt( sum(u_error.^2)/Nu )
%     v_error_2(jj,j) = sqrt( sum(v_error.^2)/Nv )
% 
%     p_error_i(jj,j) = max(abs(p_error))
%     p_error_2(jj,j) = sqrt( sum(p_error.^2)/Np)    
%     
% else

% exact solution:
% t_factor1 = exp(-4*pi^2*t_end/Re);
% t_factor2 = exp(-2*pi^2*t_end/Re);
% u_exact = -sin(pi*xu).*cos(pi*yu)*t_factor2;
% v_exact = cos(pi*xv).*sin(pi*yv)*t_factor2;
% p_exact = 0.25*(cos(2*pi*xpp)+cos(2*pi*ypp))*t_factor1;

% solution at small time step: only spatial error
sol_exact = load('../../results/TG_temporal/TG_Re100_N20_dt1e-3_Dir.mat');
V_exact = sol_exact.V;
p_exact = sol_exact.p;
u_exact = V_exact(1:Nu);
v_exact = V_exact(Nu+1:end);

u_error = uh-u_exact(:);
v_error = vh-v_exact(:);
u_error_i(j) = max(abs(u_error));
v_error_i(j) = max(abs(v_error));
u_error_2(j) = sqrt( sum(u_error.^2)/Nu );
v_error_2(j) = sqrt( sum(v_error.^2)/Nv );

% k_exact = t_factor1;
% k_error_i(j) = abs(k_exact(end)-k(end));


% Gpx_exact = -0.5*pi*sin(2*pi*xu)*t_factor1;
% Gpy_exact = -0.5*pi*sin(2*pi*yv)*t_factor1;
% Gpx_num   = Omu_inv.*(Gx*p+y_px);
% Gpy_num   = Omv_inv.*(Gy*p+y_py);
% 
% diff_px   = Gpx_num-Gpx_exact(:);
% diff_py   = Gpy_num-Gpy_exact(:);
% 
% Gp_error      = [diff_px; diff_py];
% Gp_error_i(j) = max( abs(Gp_error) )
% Gp_error_2(j) = sqrt( sum(Gp_error.^2)/(Nu+Nv) )

% figure
% contourf(xu,yu,reshape(diff_px,Nux_in,Nuy_in));
% figure
% contourf(xv,yv,reshape(diff_py,Nvx_in,Nvy_in));


p       = p - mean(p(:));
p_error      = p(:)-p_exact(:);
p_error_i(j) = max(abs(p_error));
p_error_2(j) = sqrt( sum(p_error.^2)/Np);


if (j==Nsim && Nsim>1)    
    
    figure    
    loglog(dt_list,u_error_i,'x-')
    hold on
    loglog(dt_list,v_error_i,'s-')
%     loglog(1./dt_list,k_error_i,'o-')
    loglog(dt_list,p_error_i,'d-');
    legend('u','v','p');
    grid on
    % convergence orders:
    disp('convergence orders:')
    get_slope(dt_list,u_error_i)
    get_slope(dt_list,v_error_i)
    get_slope(dt_list,p_error_i)
end

% % convection
% cu     = uh;
% cv     = vh;
% convection;
% 
% % diffusion
% d2u    = Diffu*uh + yDiffu;
% d2v    = Diffv*vh + yDiffv;
% 
% 
% du2dx_ex = 2*pi*sin(pi*xu).*cos(pi*xu).*cos(pi*yu).^2*t_factor1;
% % du2dx_error(j) = max(abs(du2dx_ex(:)-Omu_inv.*du2dx))
% duvdy_ex = -pi*cos(pi*xu).*sin(pi*xu).*(2*cos(pi*yu).^2-1)*t_factor1;
% % duvdy_error(j) = max(abs(duvdy_ex(:)-Omu_inv.*duvdy))
% convu_error(j) = max(abs(du2dx_ex(:) + duvdy_ex(:) - Omu_inv.*convu))
% 
% duvdx_ex = -pi*cos(pi*yv).*sin(pi*yv).*(2*cos(pi*xv).^2-1)*t_factor1;
% % duvdx_error(j) = max(abs(duvdx_ex(:)-Omv_inv.*duvdx))
% dv2dy_ex = 2*pi*sin(pi*yv).*cos(pi*yv).*cos(pi*xv).^2*t_factor1;
% % dv2dy_error(j) = max(abs(dv2dy_ex(:)-Omv_inv.*dv2dy))
% convv_error(j) = max(abs(duvdx_ex(:) + dv2dy_ex(:) - Omv_inv.*convv))
% 
% Diff_u = Omu_inv.*d2u;
% Diff_u_ex = 2*(pi)^2*nu*sin(pi*xu).*cos(pi*yu)*t_factor2;
% diffu_error(j) = max(abs(Diff_u(:)-Diff_u_ex(:)))
% 
% Diff_v = Omv_inv.*d2v;
% Diff_v_ex = -2*(pi)^2*nu*cos(pi*xv).*sin(pi*yv)*t_factor2;
% diffv_error(j) = max(abs(Diff_v(:)-Diff_v_ex(:)))


% figure
% surf(xp,yp,reshape(p_exact(:)-p,Npx,Npy)');
% xlabel('x')
% ylabel('y')

% keyboard
% addpath('results/Taylor_Green/');
% exact_rhs;

% end

% if (j==length(dt_list))
%     
%     figure(1)
%     loglog(dt_list,u_error_i,'bx-');
%     
%     figure(2)
%     loglog(dt_list,p_error_i,'rx-');
%     
% end