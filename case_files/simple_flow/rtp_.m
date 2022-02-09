% disp('rtp.m is empty')

%% nice velocity plot
% streamfunction
psi = get_streamfunction(V,t,options);

figure(23)
[up,vp,qp] = get_velocity(V,t,options);

subplot(2,1,1)
contour(xp,yp,up')
hold on
colorbar
% plot([2,2],[-0.5,0.5],'k-','LineWidth',3) % actuator disk
%%
max_vis = max(up,[],'all');
min_vis = min(up,[],'all');
max_psi = max(psi,[],'all');
min_psi = min(psi,[],'all');

psi_vis = min_vis + (max_vis-min_vis)*(psi-min_psi)/(max_psi-min_psi);
contour(x(2:end-1),y(2:end-1),reshape(psi_vis,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);
%%
title('u velocity component')
axis equal
hold off

subplot(2,1,2)
contour(xp,yp,vp')
hold on
colorbar
% plot([2,2],[-0.5,0.5],'k-','LineWidth',3) % actuator disk
%%
max_vis = max(vp,[],'all');
min_vis = min(vp,[],'all');
max_psi = max(psi,[],'all');
min_psi = min(psi,[],'all');

psi_vis = min_vis + (max_vis-min_vis)*(psi-min_psi)/(max_psi-min_psi);
contour(x(2:end-1),y(2:end-1),reshape(psi_vis,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);
%%
title('v velocity component')
axis equal
hold off

%% kinetic energy (actually mostly computations)


