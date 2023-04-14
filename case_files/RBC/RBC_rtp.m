% show_figures = 0;
show_figures = 1;

%% real-time plotting RBC and computing the Nusselt number

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

Nu  = options.grid.Nu;
Nv  = options.grid.Nv;
Npx = options.grid.Npx;
Npy = options.grid.Npy;

uh   = V(options.grid.indu);
vh   = V(options.grid.indv);
pres = reshape(p,Npx,Npy);
Temp = reshape(T,Npx,Npy);

[up,vp,qp] = get_velocity(V,t,options);

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

% width of lower plate
Lx   = options.grid.x2 - options.grid.x1;

% determine dTdy based on the difference operators:
% this will be first order accurate but consistent with the global balances
dTdy = options.discretization.STy*T+options.discretization.ySTy;
dTdyL = dTdy(1:Npx);
dTdyU = dTdy(end-Npx+1:end);

NusseltL = sum(-dTdyL.*hx)/Lx; % integrate over lower plate
NusseltU = sum(-dTdyU.*hx)/Lx; % integrate over upper plate
Rey_comp = sqrt(Ra*(sum(uh.*uh)/Nu+sum(vh.*vh)/Nv)/Pr); %Calculation of Reynolds number

% old implementation based on 'manually' getting derivatives 
% % get temperature derivative at lower plate
% TLo  = TBC(xp,y(1),t,options);
% % T at first grid points
% T1   = T(1:Npx);
% % T at second row of grid points
% T2   = T(Npx+1:2*Npx);
% 
% 
% % approximate derivative at lower plate with first order stencil
% dTdyL_1 = (T1 - TLo)/(0.5*hy(1));
% % approximate derivative at lower plate with second order stencil
% % assuming a uniform grid in y-dir
% dTdyL_2 = (-(1/3)*T2 + 3*T1 - (8/3)*TLo)/(hy(1));
% 
% % width of lower plate
% Lx   = options.grid.x2 - options.grid.x1;
% 
% NusseltL_1 = sum(-dTdyL_1.*hx)/Lx; % integrate over lower plate
% NusseltL_2 = sum(-dTdyL_2.*hx)/Lx; % integrate over lower plate
% 
% 
% % get temperature derivative at upper plate
% TUp  = TBC(xp,y(end),t,options);
% % T at last grid points
% T1   = T(end-Npx+1:end);
% % T at second row of grid points
% T2   = T(end-2*Npx+1:end-Npx);
% 
% 
% % approximate derivative at lower plate with first order stencil
% dTdyU_1 = -(T1 - TUp)/(0.5*hy(end));
% % approximate derivative at lower plate with second order stencil
% % assuming a uniform grid in y-dir
% dTdyU_2 = -(-(1/3)*T2 + 3*T1 - (8/3)*TUp)/(hy(end));
% 
% NusseltU_1 = sum(-dTdyU_1.*hx)/Lx; % integrate over upper plate
% NusseltU_2 = sum(-dTdyU_2.*hx)/Lx; % integrate over upper plate


% currently, store the second order approximation
NusseltL_time(n,1) = NusseltL;
NusseltU_time(n,1) = NusseltU;
Rey_time(n,1) = Rey_comp;
if show_figures
    figure(2)
    plot(t,NusseltL,'ks');
    grid on
    hold on
    plot(t,NusseltU,'kx');
    ylim([0 30])
    xlabel('t')
    ylabel('Nu')
    set(gcf,'color','w');
    set(gca,'LineWidth',1,'FontSize',14);
end
if show_figures
    figure(203)
    plot(t,Rey_comp,'ks');
    grid on
    hold on
%     ylim([0 1000])
    xlabel('t')
    ylabel('Re')
    set(gcf,'color','w');
    set(gca,'LineWidth',1,'FontSize',16);
end
%% check energy conservation properties      krishan
% change in total energy should be due to int T*v dOmega, in case viscous
% dissipation is included and no boundary contributions are present

% de_pot = vh'*(options.discretization.AT_v*T);
% diffT   = diffusion_temperature(T,t,options,0);
% de_cond = sum(diffT);
% 
% if show_figures
%     figure(3)
%     cmap = get(gca,'ColorOrder');
%     if (n>1)
%         % note that k includes e_int and (1/2)*u^2
%         plot(t,(k(n)-k(n-1))/dt,'s','Color',cmap(1,:));
%         hold on
%         plot(t,de_pot + de_cond,'o','Color',cmap(2,:));
%         grid on
%         title('d/dt (e_k + e_{int}) vs. potential energy and conduction');
% 
%         xlabel('t')
%         ylabel('energy change');
%         set(gcf,'color','w');
%         set(gca,'LineWidth',1,'FontSize',14);
%         legend('d/dt (e_k + e_{int})','potential energy source + conduction');
%     end
% end

%% create 2D plots

%% velocity
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

%% temperature    krishan
if show_figures
    figure(1)
    set(gcf,'color','w');
    % l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
    l=linspace(0,1,20);
    % l = 20;
    contour(xp,yp,Temp',l,'LineWidth',2);
    hold on
    quiver(xp,yp,up',vp');
    axis equal
    axis([x1 x2 y1 y2]);
    xlabel('x');
    ylabel('y');
    grid
    title('temperature');
    colorbar
    set(gca,'LineWidth',1,'FontSize',14);
    hold off
end
% vis_p(T,options);
% figure(501)
% plot(T);
% hold on
% xlabel('number of cells');
% ylabel('Temperature');
% disp("Hi");
%% streamfunction

% figure
% labels = 10;
% contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
%
% axis equal
% axis([x1 x2 y1 y2]);
%
% xlabel('x');
% ylabel('y');

