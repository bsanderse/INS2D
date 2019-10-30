%% real-time plotting for Taylor-Green

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% vorticity

figure(1)
omega = get_vorticity(V,t,options);
omega = reshape(omega,Nx-1,Ny-1);
labels= linspace(-6,6,20);
% contour(x,y,omega',labels);
contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);
axis square
colorbar
grid

%% velocity

% figure(1)
% [up,vp,qp] = get_velocity(V,t,options);
% qp = reshape(qp,Npx,Npy);
% labels= linspace(0,1,20);
% % labels=20;
% contour(xp,yp,qp',labels);
% axis square
% colorbar
% grid