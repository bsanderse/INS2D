%% post-processing RBC results

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    
%comment added
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
        alfa4*(NusseltU - NusseltL) - alfa3_Phi_tot
        
        % alternatively, we can compare Phi_tot directly to
        % sum(FTemp):
        FTemp = conv_diff_temperature(T,V,t,options,0);
        sum(FTemp) + alfa3_Phi_tot % 1^T * (-conv +diff)
        
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


% plot Nusselt over time
time = 0:dt:t_end;
figure
plot(time(1:rtp.n:end),NusseltL_time(1:rtp.n:end),'r')
%plot(time(1:2:end),NusseltL_time(1:2:end),'s')
grid on
hold on
plot(time(1:rtp.n:end),NusseltU_time(1:rtp.n:end),'b')
%plot(time(1:2:end),NusseltU_time(1:2:end),'x')
ylim([0 50*NusseltL])
xlabel('t')
ylabel('Nu');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);


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
