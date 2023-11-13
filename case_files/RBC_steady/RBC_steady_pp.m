%% post-processing RBC results

% line  = {'r-','b-','k-','m-','g-'};
% color = char(line(j));

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

alfa1 = options.temp.alfa1;
alfa2 = options.temp.alfa2;
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
        %% check internal energy balance
        % difference between the upper and lower Nusselt number,
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
        
        %% check thermal dissipation balance (T^2 balance)
        % note: diffusion_temperature already includes alfa4
        diffT   = diffusion_temperature(T,t,options,0);
        % note we should have T'*convT=0
        % we should also have T'*diffT = alfa*NusseltL - epsilonT
        % check this: (note that thermal_dissipation already includes
        % alfa4)
        epsilonT = thermal_dissipation(T,t,options);
        T'*diffT/Lx - (alfa4*NusseltL - epsilonT/Lx)
        % then we also have the balance
        alfa4*NusseltL - epsilonT/Lx + T'*alfa3_Phi/Lx
        % store for plotting
        epsilonT_array(j,1) = epsilonT/Lx
        TPhi_array(j,1) = T'*alfa3_Phi/Lx
        alfa4_array(j,1) = alfa4
        
        %% check kinetic energy balance
        % we need integral of Phi:
        % convert vector to 2D field to ease integration
        Phi2D = reshape(Phi,options.grid.Npx,options.grid.Npy);
        % integrate in x-direction
        Phi2D_int = sum(Phi2D,1);
        % we now have a vector which is only a function of y, which we
        % integrate "at each point", i.e. int_{0}^{y}
        Phi2D_int_cum = cumsum(Phi2D_int);
        Phi2D_int_cum_tot = sum(Phi2D_int_cum);
        % note that Phi includes alfa1 as scaling
        alfa2*alfa4*(NusseltL - 1) - sum(Phi) + (alfa2/gamma)*Phi2D_int_cum_tot*options.grid.hy(1)
end

NusseltL_array(j,1) = NusseltL
NusseltU_array(j,1) = NusseltU

% plot Nusselt over time
if (options.case.steady == 0)
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
% hold on
% quiver(xp,yp,up',vp');
axis equal
axis([x1 x2 y1 y2]);
xlabel('x');
ylabel('y');
grid
title('temperature');
colorbar
set(gca,'LineWidth',1,'FontSize',14);

%% 
if (j==Nsim && Nsim>1)
    figure
    loglog(Ra_list,NusseltL_array)
    hold on
    loglog(Ra_list,NusseltU_array)
    loglog(Ra_list,epsilonT_array./alfa4_array)
    xlim([1e3 5e5])
    ylim([0.9 8])
end

