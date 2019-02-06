%% real-time plotting lid-driven cavity 

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

%load Botella-Peyret data
% run('results/LDC/BP.m');

Nu = options.grid.Nu;
Nv = options.grid.Nv;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
pres = reshape(p,Npx,Npy);

% shift pressure to get zero pressure in the centre
if (floor(Nx/2)==Nx/2 && floor(Ny/2)==Ny/2)
    pres_ = pres-(pres(Nx/2+1,Ny/2+1)+pres(Nx/2,Ny/2))/2;
else
    pres_ = pres-pres(ceil(Nx/2),ceil(Ny/2));
end

% vorticity
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Nx-1,Ny-1);

% streamfunction
% psi = get_streamfunction(V,t,options);


%% create 2D plots

%% velocity
[up,vp,qp] = get_velocity(V,t,options);
list = linspace(0,1,20);
% list = 20;
figure(1)
set(gcf,'color','w');
pcolor(xp,yp,qp')
shading interp
% [~,c]=contour(xp,yp,qp',list);
% c.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar
caxis([0 1])
% grid
% title('velocity')
% set(gca,'LineWidth',1);


%% vorticity
% figure(1)
% labels = [-5 -4 -3 -2 -1 -0.5 0 0.5 1 2 3]; % suitable for Re=1000
% % labels = 30;
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
% figure
% 
% l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
% contour(xp,yp,pres_',l,'LineWidth',2);
% axis equal
% axis([x1 x2 y1 y2]);
% xlabeltex('x',14);
% ylabeltex('y',14);
% grid
% title('pressure');
% colorbar
% set(gca,'LineWidth',1)


%% streamfunction

% figure
% % Re=1000:
% labels = [1.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5 1e-6 0 -1e-10 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
% % Re=5000:
% % labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
% % Re=10000:
% % labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
% 
% contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
% 
% axis equal
% axis([x1 x2 y1 y2]);
% 
% xlabel('x');
% ylabel('y');

