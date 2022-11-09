%% post-processing RBC results

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

%load Botella-Peyret data
% run('results/LDC/BP.m');

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    

uh   = V(options.grid.indu);
vh   = V(options.grid.indv);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
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
disp(['residual of momentum and energy equations ' num2str(max(abs(Fres)))])

%% compute Nusselt number with final solution
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

Nu_1 = sum(-dTdy_1.*hx) % integrate over lower plate
Nu_2 = sum(-dTdy_2.*hx) % integrate over lower plate

% plot Nusselt over time
time = 0:dt:t_end;
figure
plot(time(1:rtp.n:end),Nusselt(1:rtp.n:end),'s')
grid on
hold on
ylim([0 5])
xlabel('t')
ylabel('Nu');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);


%% create 2D plots

%% vorticity
figure
% labels = [-5 -4 -3 -2 -1 -0.5 0 0.5 1 2 3]; % suitable for Re=1000
labels = 30;
contour(x(2:end-1),y(2:end-1),omega',labels,'LineWidth',2);
% 
axis equal
axis([x1 x2 y1 y2]);
% 
xlabeltex('x',14);
ylabeltex('y',14);
colorbar
grid
title('vorticity')
set(gca,'LineWidth',1);

%% pressure
figure

l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
contour(xp,yp,pres_',l,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',1)


%% streamfunction

figure
labels = 20;
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);

axis equal
axis([x1 x2 y1 y2]);

xlabel('x');
ylabel('y');


%% temperature
figure
set(gcf,'color','w');
% l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
l=linspace(-0.5,0.5,20);
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


%% convergence graphs

if (j==Nsim && Nsim>1)    
    
    figure
    if (Re==1000)
        umax_exact = 0.3885698;
        vmin_exact = -0.5270771;
        vmax_exact = 0.3769447;
    end
    if (Re==100)
        umax_exact = 0.2140424;
        vmin_exact = -0.2538030;
        vmax_exact = 0.1795728;
    end
    error_u = abs(umax_list - umax_exact);
    error_vmin = abs(vmin_list - vmin_exact);
    error_vmax = abs(vmax_list - vmax_exact);
    
    loglog(1./mesh_list,error_u,'kx-')
    hold on
    loglog(1./mesh_list,error_vmin,'ks-')
    loglog(1./mesh_list,error_vmax,'ko-')
end