%% real-time plotting for Taylor-Green

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% vorticity

% figure(1)
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Npx+1,Npy+1);
% labels= linspace(-6,6,20);
% contour(x,y,omega',labels,'LineWidth',2);
% axis square
% colorbar
% caxis([min(labels) max(labels)]);
% grid

%% velocity

figure(1)
set(gcf,'color','w');
[up,vp,qp] = get_velocity(V,t,options);
qp = reshape(qp,Npx,Npy);
labels= linspace(0,1,20);
% labels=20;
contour(xp,yp,qp',labels,'LineWidth',2);
hold on
x_skip = floor(Npx/20);
vec_x_skip = 1:x_skip:Npx;
y_skip = floor(Npy/20);
vec_y_skip = 1:y_skip:Npy;

quiver(xp(vec_x_skip),yp(vec_y_skip),up(vec_x_skip,vec_y_skip),vp(vec_x_skip,vec_y_skip));
hold off
axis([x1 x2 y1 y2]);
axis square
colorbar
caxis([min(labels) max(labels)]);
grid