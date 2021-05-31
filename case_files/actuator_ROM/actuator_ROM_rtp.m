%% real-time plotting unsteady actuator

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
% [up,vp,qp] = get_velocity(V,t,options);
% list = linspace(0.5,1.15,20);
% % list = 20;
% figure(1)
% set(gcf,'color','w');
% % pcolor(xp,yp,qp')
% [~,c]=contour(xp,yp,qp',list);
% c.LineWidth = 2;
% axis equal
% axis([x1 x2 y1 y2]);
% colorbar


% caxis([0 1])
% grid
% title('velocity')
% set(gca,'LineWidth',1);

%% ROM-FOM velocity error
% figure(1)
% clf
% Axes = zeros(2);
% [upROM,vpROM,qpROM] = get_velocity(V,t,options);
% [upFOM,vpFOM,qpFOM] = get_velocity([snapshots.uh_total(n,:) snapshots.vh_total(n,:)]',t,options);
% 
% qerror = qpROM-qpFOM;
% s = pcolor(qerror');
% % caxis manual;          % allow subsequent plots to use the same color limits
% caxis([0 0.05]);
% % axis equal
% % axis([x1 x2 y1 y2]);
% s.EdgeColor = 'none';
% title('ROM - FOM velocity error')
% % colorbar('Location','east')
% colorbar
% hold off

%% ROM-FOM velocity error componentwise
V_FOM = [snapshots.uh_total(n,:) snapshots.vh_total(n,:)]';
[upFOM,vpFOM,qpFOM] = get_velocity(V_FOM,t,options);
[upROM,vpROM,qpFOM] = get_velocity(V,t,options);

figure(1)
title('spatial error')
subplot(2,2,1)
uerror = (upFOM-upROM)';
s = pcolor(uerror);
caxis manual;          % allow subsequent plots to use the same color limits
caxis([0 0.05]);
% axis equal
% axis([x1 x2 y1 y2]);
s.EdgeColor = 'none';
title('ROM - FOM  u velocity error')
% colorbar('Location','east')
colorbar
hold off
 
% figure
subplot(2,2,3)
verror = (vpFOM-vpROM)';
s = pcolor(verror);
caxis manual;          % allow subsequent plots to use the same color limits
caxis([0 0.05]);
% axis equal
% axis([x1 x2 y1 y2]);
s.EdgeColor = 'none';
title('ROM - FOM  v velocity error')
% colorbar('Location','east')
colorbar
hold off

if t == 0
    umaxerror = uerror;
    vmaxerror = verror;
end
umaxerror = bsxfun(@max,umaxerror,uerror);
vmaxerror = bsxfun(@max,vmaxerror,verror);

subplot(2,2,2)
s = pcolor(umaxerror);
caxis manual;          % allow subsequent plots to use the same color limits
caxis([0 0.05]);
% axis equal
% axis([x1 x2 y1 y2]);
s.EdgeColor = 'none';
title('ROM - FOM  max u velocity error so far')
% colorbar('Location','east')
colorbar
hold off
 
% figure
subplot(2,2,4)
s = pcolor(vmaxerror);
caxis manual;          % allow subsequent plots to use the same color limits
caxis([0 0.05]);
% axis equal
% axis([x1 x2 y1 y2]);
s.EdgeColor = 'none';
title('ROM - FOM  max v velocity error so far')
% colorbar('Location','east')
colorbar
hold off


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
% figure
% l = 20; 
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
% labels = 20;
% contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
% 
% axis equal
% axis([x1 x2 y1 y2]);
% 
% xlabel('x');
% ylabel('y');

