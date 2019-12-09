% post-processing backward facing step

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Nu    = options.grid.Nu;
Nv    = options.grid.Nv;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;
yu     = options.grid.yu;

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);



%% velocity
[up,vp,qp] = get_velocity(V,t,options);

% BFS:
% l = [0.05 0.1 0.15 0.2 0.4 0.6 0.8 1 1.2 1.4];
list = 0.1:0.05:1.3;
% list = 25;

figure(1)
contour(xp,yp,qp',list)
axis([x1 x2 y1 y2]);
axis equal
colorbar
set(gca,'LineWidth',2)
grid
title('velocity');
colorbar
hold off

% %% pressure
% pres = reshape(p,Npx,Npy);
% 
% %BFS:
% % set the pressure to zero at the corner of the step
% pres  = pres-interp2(xp,yp,pres',0,0,'spline');
% shift = abs((pres(1,Npy/2)+pres(1,Npy/2+1))/2);
% % shift = abs(pres(1,Npy/2));
% pres = pres + shift;
% 
% % BFS:
% l = [0.01:0.01:0.1 0.12:0.02:0.24];
% 
% figure
% contour(xp,yp,pres',l,'LineWidth',2);
% axis equal
% axis([x1 x2 y1 y2]);
% xlabeltex('x',14);
% ylabeltex('y',14);
% grid
% title('pressure');
% colorbar
% set(gca,'LineWidth',2)
% 
% %% vorticity
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Nx-1,Ny-1);
% 
% figure
% labels = -8:2:10;
% contour(x(2:end-1),y(2:end-1),reshape(omega,Nx-1,Ny-1)',labels,'LineWidth',1);
% % 
% axis equal
% axis([x1 x2 y1 y2]);
% % 
% xlabeltex('x',14);
% ylabeltex('y',14);
% colorbar
% grid
% title('vorticity')
% set(gca,'LineWidth',2);
% % hold off;
% 
% %% streamfunction
% psi = get_streamfunction(V,t,options);
% labels = 25;
% contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
% axis equal
% axis([x1 x2 y1 y2]);
% xlabel('x');
% ylabel('y');
% colorbar
% grid
% title('streamfunction')
% set(gca,'LineWidth',2)

