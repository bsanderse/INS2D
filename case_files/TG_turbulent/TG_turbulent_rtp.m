%% real-time plotting for Taylor-Green

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% vorticity

figure(1)
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx+1,Npy+1);
labels= linspace(-6,6,20);
contour(x,y,omega',labels,'LineWidth',2);
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