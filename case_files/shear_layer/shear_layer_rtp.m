%% real-time plotting for shear layer

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% vorticity
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
labels= linspace(-2,2,20);
[~,c] = contour(x,y,omega',labels);
c.LineWidth = 1;
axis square
colorbar
% grid