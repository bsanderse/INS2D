%% real-time plotting for shear layer

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% vorticity
% % compare e.g. with PhD thesis figure 3.4
% % set(0,'defaultlinelinewidth',3)
% 
% figure(1)
% set(gcf,'color','w');
% % set(gca,'LineWidth',6);
% ax = gca;
% ax.LineWidth = 2;
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Npx+1,Npy+1);
% % for Re=1000: labels = -4:0.5:4;
% labels= linspace(-2,2,20);
% [~,c] = contour(x,y,omega',labels);
% c.LineWidth = 1;
% axis square
% colorbar
% title('vorticity')
% xlabel('x dimension')
% ylabel('y dimension')
% % grid

%% velocity
Nu  = options.grid.Nu;
Nv  = options.grid.Nv;
Npx = options.grid.Npx;
Npy = options.grid.Npy;   

% streamfunction
% psi = get_streamfunction2(V,t,options);

% uh   = V(1:Nu);
% vh   = V(Nu+1:Nu+Nv);
% pres = reshape(p,Npx,Npy);
% 
% figure(1)
% clf
% Axes = zeros(2);
% [up,vp,qp] = get_velocity(V,t,options);
% % list = 20;
% figure(1)
% set(gcf,'color','w');
% 
% % Axes(1) = axes;
% % contour(Axes(1),x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);
% % % colorbar
% % axis equal
% % axis([x1 x2 y1 y2]);
% % hold off
% 
% Axes(2) = axes;
% % pcolor(xp,yp,qp')
% list = linspace(0.6,1.1,20);
% [~,c]=contour(Axes(2),xp,yp,qp');%,list);
% c.LineWidth = 1;
% axis equal
% axis([x1 x2 y1 y2]);
% colorbar('Location','east')
% title('velocity magnitude')


% linkaxes(Axes)
hold off

figure(2)
clf
figure(2)
set(gcf,'color','w');
[up,vp,qp] = get_velocity(V,t,options);

subplot(1,2,1)
[~,c_u]=contour(xp,yp,up');
c_u.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar('Location','east')
title('velocity u component')


subplot(1,2,2)
[~,c_v]=contour(xp,yp,vp');
c_v.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar('Location','east')
title('velocity v component')

hold off