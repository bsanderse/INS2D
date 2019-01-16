function [up,vp,qp] = get_velocity(V,t,options)
% velocity values at pressure points

% restrict to velocities at pressure positions, using BMx and BMy,
% constructed in operator_divergence

% evaluate boundary conditions at current time
% boundary_conditions;
% interpolate_bc;
% operator_bc_divergence;
% operator_bc_momentum;
options = set_bc_vectors(t,options);

Au_ux  = options.discretization.Au_ux;
yAu_ux = options.discretization.yAu_ux;
Av_vy  = options.discretization.Av_vy;
yAv_vy = options.discretization.yAv_vy;
Bup = options.discretization.Bup;
Bvp = options.discretization.Bvp;
Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nu  = options.grid.Nu;
Nv  = options.grid.Nv;

uh = V(1:Nu);
vh = V(Nu+1:Nu+Nv);

uh_temp = uh; vh_temp = vh;
if (options.case.ibm==1)
    % blank inside values for plotting purposes, in case of immersed boundary method
    % the real values should NOT be blanked, they are needed in the
    % ghost-cell method.
    uh_temp(inside_incl_interface_u(:)==1)=0;
    vh_temp(inside_incl_interface_v(:)==1)=0;
end

up    = reshape( Bup*(Au_ux * uh_temp + yAu_ux), Npx, Npy);
vp    = reshape( Bvp*(Av_vy * vh_temp + yAv_vy), Npx, Npy);

qp    = sqrt(up.^2+vp.^2); 


% with these velocities we can make plots like:
% contour(xp,yp,qp')


 
%% BFS:
% % l = [0.05 0.1 0.15 0.2 0.4 0.6 0.8 1 1.2 1.4];
% % list = 0.1:0.05:1.3;
% list = 25;
% 
% figure(10)
% % pcolor(xp,yp,qp')
% contour(xp,yp,qp',list)
% % surf(xp,yp,qp');
% % contour(xp,yp,up',25);
% % shading interp
% axis([x1 x2 y1 y2]);
%  axis equal
% % axis tight
% % axis square
% colorbar
% 
% % print(gcf,'-depsc',['results/nonaligned_convection/surf_' num2str(Nx)]);
% 
% xlabel('$x$','Interpreter','Latex');
% ylabel('$y$','Interpreter','Latex');


%% get wake profiles
% u = reshape(uh,Nux_in,Nuy_in);
% v = reshape(vh,Nvx_in,Nvy_in);
% uwake1 = interp2(xin',yp,u',x_c+0.5,yp);
% uwake2 = interp2(xin',yp,u',x_c+5,yp);
% vwake1 = interp2(xp',yin,v',x_c+0.5,yin);
% vwake2 = interp2(xp',yin,v',x_c+5,yin);
% 
% 
% figure(2)
% plot(uwake1,yp,'rx-')
% hold on
% plot(uwake2,yp,'ro-')
% figure(3)
% plot(vwake1,yin,'bx-')
% hold on
% plot(vwake2,yin,'bo-')


% figure(2)
% contour(xp,yp,qp',25)
% axis equal
% axis tight
% axis([x1 x2 y1 y2])
% colorbar
% print(gcf,'-depsc',['results/nonaligned_convection/contour_' num2str(Nx)]);


% actuator disk
% hold on
% h1=plot([x(xmid+1) x(xmid+1)],[yp(yrange(1)) yp(yrange(end))],'k-');
% set(h1,'LineWidth',1.5);
% h2=plot([x(2*xmid+1) x(2*xmid+1)],[yp(yrange(1)+floor(nd/2)) yp(yrange(end)+floor(nd/2))],'k-');
% set(h2,'LineWidth',1.5);

% airfoil
% hold on
% plot(x_k,y_k,'bx-');

% hold off;
