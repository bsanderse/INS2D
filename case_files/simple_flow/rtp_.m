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
if options.verbosity.energy_verbosity == 1
    NV = options.grid.NV;
    Np = options.grid.Np;
    
    Q_h = options.discretization.Q_h;
    visc = options.case.visc;
    
    k_diff = - visc*norm(Q_h*V);
    
    K_h = options.discretization.K_h;
    I_h = options.discretization.I_h;
    A_h = options.discretization.A_h;
    y_A = options.discretization.y_A;
    y_I = options.discretization.y_I;
    
    k_conv = -V'*K_h*(spdiags(I_h*V+y_I,0,NV,NV)*(A_h*V+y_A));
    
    p_h = zeros(Np,1);
    p_h = pressure_additional_solve(V,p_h,t,options);
    
    M_h = options.discretization.M;
    
    k_pres = (M_h*V)'*p_h;
    
    y_G = [options.discretization.y_px; ...
        options.discretization.y_py];
    
    y_D = [options.discretization.yDiffu; ...
        options.discretization.yDiffv];

    k_presBC = 
    
    k_diffBC = 
end