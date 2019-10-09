%% real-time plotting for shear layer

Npx = options.grid.Npx;
Npy = options.grid.Npy;

yu = options.grid.yu;
yv = options.grid.yv;

d     = options.fluid.d_layer;
U1    = options.fluid.U1;
U2    = options.fluid.U2;


u_base = 0.5*(U1 + U2) + 0.5*(U1-U2)*tanh(yu./d);
v_base = 0*yv;
V_pert = V - [u_base(:);v_base(:)];
% mode_u = options.fluid.mode_u;
% omega  = options.fluid.omega;
% alpha  = options.fluid.alpha;    
% pert   = options.fluid.pert;
% 
% I      = sqrt(-1);
% u_pert = real(mode_u(yu).*exp(I*(alpha*xu - omega*t)));
% v_pert = real(mode_v(yv).*exp(I*(alpha*xv - omega*t)));
% 
% u = u + pert*u_pert;
% v = v + pert*v_pert;
%%
% [up,vp,qp] = get_velocity(V,t,options);
% list = linspace(0.7,1.15,20);
% % list = 20;
% figure(1)
% set(gcf,'color','w');
% % pcolor(xp,yp,qp')
% [~,c]=contour(xp,yp,qp',list);
% c.LineWidth = 2;
% axis equal
% axis([x1 x2 y1 y2]);
% colorbar
% % caxis([0 1])
% % grid
% % title('velocity')
% % set(gca,'LineWidth',1);

%% vorticity
set(0,'defaultlinelinewidth',3)

figure(1)
set(gcf,'color','w');
% set(gca,'LineWidth',6);
ax = gca;
ax.LineWidth = 2;
omega = get_vorticity(V_pert,t,options);

omega = reshape(omega,Npx-1,Npy-1);
% for Re=1000: labels = -4:0.5:4;
% labels= linspace(-0.5,0.2,20);
labels=20;
[~,c] = contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);

% [~,c] = contour(x,y,omega',labels);
c.LineWidth = 1;
axis square
colorbar
% grid

%% pressure
% figure(1)
% set(gcf,'color','w');
% 
% pres = reshape(p,Npx,Npy);
% 
% l = 20; 
% % l = linspace(-5e-3,5e-3,20);
% contour(xp,yp,pres',l,'LineWidth',2);
% axis equal
% axis([x1 x2 y1 y2]);
% xlabeltex('x',14);
% ylabeltex('y',14);
% grid
% title('pressure');
% colorbar
% % caxis([l(1) l(end)]);
% set(gca,'LineWidth',1)