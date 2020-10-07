%% post-processing actuator disk of Maarten vd Broek's testcase


%% get grid variables
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Nu     = options.grid.Nu;
Nv     = options.grid.Nv;
Npx = options.grid.Npx;
Npy = options.grid.Npy;
xin  = options.grid.xin;
yin  = options.grid.yin;



%% get solution variables
uh   = V(1:Nu);
vh   = V(Nu+1:end);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);

figure
contour(xp,yp,pres');

%% velocity
[up,vp,qp] = get_velocity(V,t,options);
% list = linspace(0,1,20);
list = 20;
figure
set(gcf,'color','w');
% pcolor(xp,yp,qp')
[~,c]=contourf(xp,yp,qp',list);
c.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar
% caxis([0 1])
% grid
% title('velocity')
% set(gca,'LineWidth',1);

%% velocity as pcolor
figure
set(gcf,'color','w');
% pcolor(xp,yp,qp')
pcolor(xp,yp,qp');
shading interp
% c.LineWidth = 1;
box on
axis equal
axis([x1 x2 y1 y2]);
colorbar

%% plot settings
% line = {'r-','b-','m-','k-','g-','y-'};
% color = char(line(j));

%% turbulent diffusion
[S11,S12,S21,S22,S_abs] = strain_tensor(V,t,options,0);

nu_t = turbulent_viscosity(S_abs,options);

% list = linspace(0,1,20);
list = 20;
figure
set(gcf,'color','w');
% pcolor(xp,yp,qp')
[~,c]=contourf(xp,yp,reshape(nu_t,Npx,Npy)',list);
c.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar