%% real-time plotting for shear layer

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% decaying_turbulence
% compare e.g. with PhD thesis figure 3.4
% set(0,'defaultlinelinewidth',3)

figure(1)
set(gcf,'color','w');
% set(gca,'LineWidth',6);
ax = gca;
ax.LineWidth = 2;
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx+1,Npy+1);
% for Re=1000: labels = -4:0.5:4;
% labels= linspace(-2,2,20);
labels = 16;
[~,c] = contour(x,y,omega',labels);
c.LineWidth = 2;
axis square
colorbar
% grid

% figure(1)
% set(gcf,'color','w');
% [up,vp,qp] = get_velocity(V,t,options);
% qp = reshape(qp,Npx,Npy);
% % labels= linspace(0,2,20);
%  labels=20;
% contour(xp,yp,qp',labels,'LineWidth',2);
% axis square
% colorbar