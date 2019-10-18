%% unsteady actuator
% see figure 6.6 in PhD thesis

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

Nu  = options.grid.Nu;
Nv  = options.grid.Nv;
Npx = options.grid.Npx;
Npy = options.grid.Npy;    

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
pres = reshape(p,Npx,Npy);

% shift pressure to get zero pressure in the centre
% if (floor(Nx/2)==Nx/2 && floor(Ny/2)==Ny/2)
%     pres_ = pres-(pres(Nx/2+1,Ny/2+1)+pres(Nx/2,Ny/2))/2;
% else
%     pres_ = pres-pres(ceil(Nx/2),ceil(Ny/2));
% end

% vorticity
omega = get_vorticity(V,t,options);
omega = reshape(omega,Nx-1,Ny-1);

% streamfunction
psi = get_streamfunction(V,t,options);


%% create 2D plots

%% velocity
figure(1)
clf
Axes = zeros(2);
[up,vp,qp] = get_velocity(V,t,options);
% list = 20;
figure(1)
set(gcf,'color','w');
Axes(1) = axes;
contour(Axes(1),x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);
axis equal
axis([x1 x2 y1 y2]);
xlabel('x/D')
ylabel('y/D')
title('velocity field')
set(gca,'FontSize',12);

hold off

Axes(2) = axes;
% pcolor(xp,yp,qp')
list = linspace(0.6,1.1,20);
[~,c]=contour(Axes(2),xp,yp,qp',list);
c.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar('Location','east')
% caxis([0 1])
% grid
% set(gca,'LineWidth',1);
% hold on
% labels=20;
hold on
plot([2,2],[-0.5,0.5],'k-','LineWidth',3) % actuator disk
hold off
set(Axes(2), 'visible', 'off');

linkaxes(Axes)

%% vorticity
% figure(1)
% labels = 30;
% contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);
% % 
% axis equal
% axis([x1 x2 y1 y2]);
% % 
% xlabeltex('x',14);
% ylabeltex('y',14);
% colorbar
% grid
% title('vorticity')
% set(gca,'LineWidth',1);

%% pressure
figure
l = 20; 
contour(xp,yp,pres',l,'LineWidth',1);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',1)


%% streamfunction

% figure
% labels = 20;
% contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
% 
% axis equal
% axis([x1 x2 y1 y2]);
% 
% xlabel('x');
% ylabel('y');

%% kinetic energy
figure
plot(0:dt:t_end,k/k(1));
title('normalised kinetic energy')

%% divergence of velocity field
figure
semilogy(0:dt:t_end,maxdiv);
title('divergence of velocity field')
