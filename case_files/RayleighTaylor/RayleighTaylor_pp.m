%% post-processing RBC results

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    

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
omega = get_vorticity(V,t,options);
omega = reshape(omega,Nx-1,Ny-1);
% streamfunction
psi = get_streamfunction(V,t,options);

%% check whether steady state is achieved

Fres = F(V,V,p,T,t,options,0,0);
disp(['residual of momentum and energy equations: ' num2str(Fres)])

%% compute Nusselt number with final solution

alfa3 = options.temp.alfa3;
alfa4 = options.temp.alfa4;

% width of lower plate
Lx   = options.grid.x2 - options.grid.x1;

% determine dTdy based on the difference operators:
% this will be first order accurate but consistent with the global balances
dTdy = options.discretization.STy*T+options.discretization.ySTy;
dTdyL = dTdy(1:Npx);
dTdyU = dTdy(end-Npx+1:end);

NusseltL = sum(-dTdyL.*hx)/Lx % integrate over lower plate
NusseltU = sum(-dTdyU.*hx)/Lx % integrate over upper plate


% alternative below is to manually get derivatives
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

% NusseltL_1 = sum(-dTdyL.*hx)/Lx % integrate over lower plate
% NusseltL_2 = sum(-dTdyL_2.*hx)/Lx % integrate over lower plate


% % get temperature derivative at upper plate
% TUp  = TBC(xp,y(end),t,options);
% % T at last grid points
% T1   = T(end-Npx+1:end);
% % T at second row of grid points
% T2   = T(end-2*Npx+1:end-Npx);
% 
% 
% % approximate derivative at upper plate with first order stencil
% dTdyU_1 = (TUp - T1)/(0.5*hy(end));
% % approximate derivative at upper plate with second order stencil
% % assuming a uniform grid in y-dir
% dTdyU_2 = ((1/3)*T2 - 3*T1 + (8/3)*TUp)/(hy(end));

% NusseltU_1 = sum(-dTdyU_1.*hx)/Lx % integrate over upper plate
% NusseltU_2 = sum(-dTdyU_2.*hx)/Lx % integrate over upper plate

switch options.temp.incl_dissipation
    case 1
        % check difference between the upper and lower Nusselt number,
        % this should equal the dissipation
        [Phi] = dissipation(V,t,options,0);
        % the computed dissipation is basically V'*D*V, which includes
        % alfa1 as scaling
        % we are interested in alfa3*||gradV||^2, so we divide by gamma
        gamma   = options.temp.gamma;
        alfa3_Phi     = (1/gamma) * Phi;
        alfa3_Phi_tot = sum(alfa3_Phi);
        alfa4*(NusseltU - NusseltL) - alfa3_Phi_tot/Lx
        
        % alternatively, we can compare Phi_tot directly to
        % sum(FTemp):
        FTemp = conv_diff_temperature(T,V,t,options,0);
        sum(FTemp)/Lx + alfa3_Phi_tot/Lx % 1^T * (-conv +diff)
        
        %% check thermal dissipation balance
        % note: diffusion_temperature already includes alfa4
        diffT   = diffusion_temperature(T,t,options,0);
        % note we should have T'*convT=0
        % we should also have T'*diffT = alfa*NusseltL - epsilonT
        % check this:
        epsilonT = thermal_dissipation(T,t,options);
        T'*diffT/Lx - (alfa4*NusseltL - epsilonT/Lx)
        % then we also have the balance
        alfa4*NusseltL - epsilonT/Lx + T'*alfa3_Phi/Lx

end

NusseltL_array(j,1) = NusseltL
NusseltU_array(j,1) = NusseltU

% plot Nusselt over time
time = 0:dt:t_end;
figure(101)
plot(time(1:rtp.n:end),NusseltL_time(1:rtp.n:end),'-')
grid on
hold on
% plot(time(1:rtp.n:end),NusseltU_time(1:rtp.n:end),'x')
ylim([0 4])
xlabel('t')
ylabel('Nu');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);


% plot energy balances
figure(102)
cmap = get(gca,'ColorOrder');

% note that k includes both e_int and (1/2)*u^2 !!
plot(time(1:rtp.n:end),dkdt(:,j),'s','Color',cmap(j,:));
hold on
plot(time(1:rtp.n:end),de_pot(:,j) + de_cond(:,j),'o','Color',cmap(j,:));

grid on
title('d/dt (e_k + e_{i}) vs. potential energy and conduction');        
xlabel('t')
ylabel('energy change');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);
legend('d/dt (e_{k} + e_{i})','potential energy source + conduction');


figure(103)
cmap = get(gca,'ColorOrder');
semilogy(time(1:rtp.n:end),abs(energy_diff(:,j)),'s','Color',cmap(j,:));
hold on

grid on
title('d/dt (e_{k} + e_{i}) - potential energy');

xlabel('t')
ylabel('energy change');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);
if (j==Nsim)
    legend('Ge=0.1, OL','Ge=1, OL','Ge=0.1, IM','Ge=1, IM');
    % ylabel('$\varepsilon_{E}$','Interpreter','Latex');
end

figure(104)
% average temperature on domain

cmap = get(gca,'ColorOrder');

plot(time(1:rtp.n:end),T_avg(:,j),'o','Color',cmap(j,:));
hold on
plot(time(1:rtp.n:end),u_avg(:,j),'x','Color',cmap(j,:));
grid on
ylabel('average temperature and average velocity');
%     ylim([0.48 0.53]);
xlim([options.time.t_start options.time.t_end]);
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);

% do some postprocessing after last run finishes
if (j==Nsim)
    legend('T, Ge=0.1, OL','u, Ge=0.1, OL','T, Ge=1, OL','u, Ge=1, OL','T, Ge=0.1, IM','u, Ge=0.1, IM','T, Ge=1, IM','u, Ge=1, IM');
    hline = findobj(gcf, 'type', 'line');
    % set(hline,'Marker','none');
    % set(hline(1:2:end),'LineStyle','--')
    % set(hline(2:2:end),'LineStyle','-')
    set(hline(5:8),'LineStyle','-')
    set(hline(1:4),'LineStyle','none','Marker','o','MarkerIndices',1:20:1001)
    set(hline(5),'Color',cmap(8,:))
    set(hline(8),'Color',cmap(1,:))
    set(hline(7),'Color',cmap(2,:))
    set(hline(6),'Color',cmap(3,:))
    set(hline(1),'Color',cmap(4,:))
    set(hline(2),'Color',cmap(3,:))
    set(hline(3),'Color',cmap(2,:))
    set(hline(4),'Color',cmap(1,:))
end

%% create 2D plots

%% vorticity
figure
set(gcf,'color','w');
labels = 20;
contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);
% 
axis equal
axis([x1 x2 y1 y2]);
% 
xlabel('x');
ylabel('y');
colorbar
grid
title('vorticity')
set(gca,'LineWidth',1);

%% pressure
figure
set(gcf,'color','w');
labels = 20; 
contour(xp,yp,pres_',labels,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabel('x');
ylabel('y');
grid
title('pressure');
colorbar
set(gca,'LineWidth',1)


%% streamfunction

figure
set(gcf,'color','w');
labels = 20;
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
title('streamfunction');
xlabel('x');
ylabel('y');
colorbar
set(gca,'LineWidth',1)

%% temperature
figure
set(gcf,'color','w');
l=linspace(0,1,20);
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
set(gca,'LineWidth',1)
