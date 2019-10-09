%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;
Nu  = options.grid.Nu;
Nv  = options.grid.Nv;

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
pres = reshape(p,Npx,Npy);

%% kinetic energy
figure
plot(time,k/k(1),'s-')
grid
title('kinetic energy')

%% vorticity
% compare e.g. with PhD thesis figure 3.4

figure
set(gcf,'color','w');

omega = get_vorticity(V,t,options);

omega = reshape(omega,Npx-1,Npy-1);
% for Re=1000: labels = -4:0.5:4;
labels= linspace(-0.5,0.2,20);
[~,c] = contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);

% axis square
colorbar
grid
title('vorticity')
set(gca,'LineWidth',1)


%%

figure
set(gcf,'color','w');

l = 20; 
contour(xp,yp,pres',l,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',1)

%% velocity
[up,vp,qp] = get_velocity(V,t,options);
% list = linspace(0.7,1.15,20);
list = 20;
figure
set(gcf,'color','w');
% pcolor(xp,yp,qp')
[~,c]=contour(xp,yp,qp',list);
c.LineWidth = 2;
axis equal
axis([x1 x2 y1 y2]);
colorbar
% caxis([0 1])
% grid
title('velocity')
set(gca,'LineWidth',1);