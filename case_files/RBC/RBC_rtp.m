%% real-time plotting RBC and computing the Nusselt number

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
Temp = reshape(T,Npx,Npy);

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

%% compute Nusselt number
% get temperature derivative at lower plate
TLo  = TBC(xp,y(1),t,options);
% T at first grid points
T1   = T(1:Npx);
% T at second row of grid points
T2   = T(Npx+1:2*Npx);


% approximate derivative at lower plate with first order stencil
dTdy_1 = (T1 - TLo)/(0.5*hy(1));
% approximate derivative at lower plate with second order stencil
% assuming a uniform grid in y-dir
dTdy_2 = (-(1/3)*T2 + 3*T1 - (8/3)*TLo)/(hy(1));

Nusselt_1 = sum(-dTdy_1.*hx); % integrate over lower plate
Nusselt_2 = sum(-dTdy_2.*hx); % integrate over lower plate

Nusselt(n,1) = Nusselt_2;

figure(2)
plot(t,Nusselt_2,'ks');
grid on
hold on
% ylim([0 5])

%% check energy conservation properties
% change in total energy should be due to int T*v dOmega
de_pot = vh'*(options.discretization.AT_v*T);
max(abs(vh));
% pause
figure(3)
if (n>1)
    % note that k includes e_int and (1/2)*u^2
    plot(t,(k(n)-k(n-1))/dt,'ks');
    hold on
    plot(t,de_pot,'bs');
    grid on
    title('dk/dt and e_{pot}');
end

%% create 2D plots

%% velocity
[up,vp,qp] = get_velocity(V,t,options);
% % list = linspace(0,1,20);
% list = 20;
% figure(1)
% set(gcf,'color','w');
% pcolor(xp,yp,qp')
% shading interp
% hold on
% quiver(xp,yp,up',vp');
% % [~,c]=contour(xp,yp,qp',list);
% % c.LineWidth = 1;
% axis equal
% axis([x1 x2 y1 y2]);
% colorbar
% caxis([0 0.2])
% hold off
% % grid
% % title('velocity')
% % set(gca,'LineWidth',1);


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

%% temperature
figure(1)
set(gcf,'color','w');
% l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
% l=linspace(-0.5,0.5,20);
l = 20;
contour(xp,yp,Temp',l,'LineWidth',2);
hold on
quiver(xp,yp,up',vp');
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('temperature');
colorbar
set(gca,'LineWidth',1)
hold off


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

