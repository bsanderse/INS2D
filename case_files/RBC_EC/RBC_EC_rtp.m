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

% show_figures = 1;
show_energy = 0;
show_Nusselt = 0;
show_temperature = 1;

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
% height (should be one)
Ly   = options.grid.y2 - options.grid.y1;

% determine dTdy based on the difference operators:
% this will be first order accurate but consistent with the global balances
dTdy = options.discretization.STy*T+options.discretization.ySTy;
dTdyL = dTdy(1:Npx);
dTdyU = dTdy(end-Npx+1:end);

NusseltL = sum(-dTdyL.*hx)/Lx; % integrate over lower plate
NusseltU = sum(-dTdyU.*hx)/Lx; % integrate over upper plate

% % old implementation based on 'manually' getting derivatives
%
% % get temperature derivative at lower plate
% TLo  = TBC(xp,y(1),t,options);
% % T at first grid points
% T1   = T(1:Npx);
% % T at second row of grid points
% T2   = T(Npx+1:2*Npx);
%
% % approximate derivative at lower plate with first order stencil
% dTdy_1 = (T1 - TLo)/(0.5*hy(1));
% % approximate derivative at lower plate with second order stencil
% % assuming a uniform grid in y-dir
% dTdy_2 = (-(1/3)*T2 + 3*T1 - (8/3)*TLo)/(hy(1));
%
% Nusselt_1 = sum(-dTdy_1.*hx)/Lx; % integrate over lower plate
% Nusselt_2 = sum(-dTdy_2.*hx)/Lx; % integrate over lower plate
%
% % we store and plot the second order approximation
% Nusselt(n,1) = Nusselt_2;

%store the Nusselt number
NusseltL_time(n,1) = NusseltL;
NusseltU_time(n,1) = NusseltU;

if (show_Nusselt)
    
    figure(2)
    plot(t,Nusselt_2,'ks');
    grid on
    hold on
    ylim([0 5])
    xlabel('t')
    ylabel('Nu')
    set(gcf,'color','w');
    set(gca,'LineWidth',1,'FontSize',14);
    
end

%% check energy conservation properties
% change in total energy should be due to int T*v dOmega, in case viscous
% dissipation is included and no boundary contributions are present

if (show_energy)
    
        figure(3)
        cmap = get(gca,'ColorOrder');
        % de_pot = vh'*(options.discretization.AT_v*T);
        % in implicit midpoint, the contribution equals
        % V^{n+1/2} * A * T^{n+1/2}
        V_mid = 0.5*(V + Vn);
        T_mid = 0.5*(T + Tn);
        % average temperature on domain
        T_avg = sum(options.grid.Omp.*T)/(Lx*Ly);
        de_pot = V_mid(options.grid.indv)'*(options.discretization.AT_v*T_mid);
        dkdt   = (k(n)-k(n-1))/dt;
        % note that k includes e_int and (1/2)*u^2

        yyaxis left
        plot(t,dkdt,'s'); %'Color',cmap(1,:));
        hold on
        plot(t,de_pot,'x');%,'Color',cmap(2,:));
        ylabel('energy change');
        ylim([-0.01 0.01]);


        yyaxis right
        plot(t,T_avg,'o'); %'Color',cmap(3,:));
        ylabel('average temperature');
        ylim([0 0.1]);

    
%         figure(4)
%         cmap = get(gca,'ColorOrder');
%         DeltaE = abs(de_pot - dkdt);
%         
%         % note that k includes e_int and (1/2)*u^2
%         semilogy(t,DeltaE,'s','Color',cmap(1,:));
%         grid on
%         hold on
%         title('de/dt - h_{pot}');
%         xlabel('t')
%         ylabel('energy change');
%         set(gcf,'color','w');
%         set(gca,'LineWidth',1,'FontSize',14);
%         %     legend('de/dt','h_{pot}');
        
%     end
    
end
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

%% temperature
if (show_temperature)

    if (n==1)

    fig1 = figure(1);
    sgtitle(['Pr = ' num2str(options.temp.Pr) ', Ra = ' num2str(options.temp.Ra) ', Ge = ' num2str(options.temp.Ge)])    
    fig1.Position(3) = [1120];
    cmap = get(gca,'ColorOrder');

    elseif (n>1)

    subplot(1,2,1)
    T_avg = sum(options.grid.Omp.*T)/(Lx*Ly);
    set(gcf,'color','w');
    plot(t,T_avg,'o','Color',cmap(3,:),'MarkerSize',10,'MarkerFaceColor',cmap(3,:));
    ylabel('average temperature');
    
%     axis equal
    axis([x1 x2 y1 y2]);
    xlabel('t');
    ylabel('T');
    xlim([0 options.time.t_end]);
    ylim([0 1])
    grid on
    title('average temperature');
    colorbar
    set(gca,'LineWidth',1,'FontSize',14);
    hold on

    subplot(1,2,2)
    V_mid = 0.5*(V + Vn);
    T_mid = 0.5*(T + Tn);
    % average temperature on domain
    de_pot = V_mid(options.grid.indv)'*(options.discretization.AT_v*T_mid);
    dkdt   = (k(n)-k(n-1))/dt;
    % note that k includes e_int and (1/2)*u^2
    
%     yyaxis left
    plot(t,dkdt,'s','Color',cmap(1,:),'MarkerSize',10);
    hold on
    plot(t,de_pot,'x','Color',cmap(2,:),'MarkerSize',10);
    ylabel('energy change');
    ylim([-0.01 0.02]);
    
%     yyaxis right
%     ylim([0 0.1]);
%     hold on

    grid on
    title('energy change'); %d/dt (e_k + e_{i}) vs. potential energy and conduction');
    xlim([0 options.time.t_end]);
    xlabel('t')
    set(gcf,'color','w');
    set(gca,'LineWidth',1,'FontSize',14);
%     if (n==1+options.visualization.rtp_n)
        legend('d/dt (e_k + e_{i})','potential energy'); %,'average temperature');
%     end
%     subplot(1,2,2)
%     plot(t,NusseltL,'rs');
%     grid on
%     hold on
%     plot(t,NusseltU,'bx');
% 	axis square
%     xlim([0 options.time.t_end]);
%     ylim([0 5])
%     xlabel('t')
%     ylabel('Nu')
%     title('Nusselt number');
%     set(gcf,'color','w');
%     set(gca,'LineWidth',1,'FontSize',14);
%     legend('Lower (hot)','Upper (cold)');
%     
%     
    end
    
end

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

