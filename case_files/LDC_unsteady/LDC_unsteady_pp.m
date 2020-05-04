%% post-processing Lid-driven cavity results

line = {'r-','b-','m-','k-','g-','y-'};
color = char(line(mod(j-1,length(line))+1));

% write output along centerlines to datafiles:
write_files = 0;
show_plots = 0;

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


%% create plots

if (show_plots == 1)
    
    %% vorticity
    figure
    labels = [-5 -4 -3 -2 -1 -0.5 0 0.5 1 2 3]; % suitable for Re=1000
    % labels = 30;
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
    % Re=1000:
    labels = [1.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5 1e-6 0 -1e-10 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
    % Re=5000:
    % labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
    % Re=10000:
    % labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
    
    contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
    
    axis equal
    axis([x1 x2 y1 y2]);
    
    xlabel('x');
    ylabel('y');
    
    %% create centerline plots
    
    if (floor(Nx/2)==Nx/2) % Nx even
        disp('min u along vertical centerline');
        
        umid = u(Nx/2,:);
        [umin, pos] = min(umid);
        disp(umin);
        disp(yp(pos));
        
        %reconstruct data with spline for more accurate result:
        if (min(hy) == max(hy)) % uniform grid
            min_y = deltay/10;
        else
            min_y = min(hy);
        end
        cs = spline(yp,umid);
        L_y   = options.grid.y2 - options.grid.y1;
        yfine = linspace(yp(1),yp(end),L_y/min_y); % linear distribution with minimum grid size
        ufine = ppval(cs,yfine);
        [umin, pos] = min(ufine);
        umax_list(j) = umin;
        ymax_list(j) = yfine(pos);
        
        %%
        figure
        hold on
        %     if (j==1)
        %     plot(uBP,yBP,'bx');
        %     end
        %     plot(-uGhia,yBP,'rx');
        plot(u(Nx/2,:),yp,color)
        xlabel('u')
        ylabel('y')
        box
        title('horizontal velocity at vertical centerline');
        
        
        %%
        figure
        hold on
        %     plot(wvBP,yBP,'bx');
        plot(-omega(Nx/2,:),y(2:end-1),color)
        xlabel('w')
        ylabel('y')
        box
        title('vorticity at vertical centerline');
        
        
        figure
        hold on
        %     if (j==1)
        %     plot(pvBP,yBP,'bx');
        %     end
        py = (pres_(Nx/2,:)+pres_(Nx/2+1,:))/2;
        plot(py,yp,color)
        xlabel('p')
        ylabel('y')
        box
        title('pressure at vertical centerline');
        
    else
        [umin pos] = min((u(ceil(Nx/2),:)+u(floor(Nx/2),:))/2);
        disp(umin);
        disp(yp(pos));
    end
    
    if (floor(Ny/2)==Ny/2) % Ny even
        disp('max and min v along horizontal centerline');
        
        vmid  = v(:,Ny/2);
        [vmax pos] = max(vmid);
        disp(vmax);
        disp(xp(pos));   %    disp(1-xp(pos));
        
        %reconstruct data with spline for more accurate results:
        if (min(hx) == max(hx)) % uniform grid
            min_x = deltax/10;
        else
            min_x = min(hx);
        end
        cs = spline(xp,vmid);
        L_x   = options.grid.x2 - options.grid.x1;
        xfine = linspace(xp(1),xp(end),L_x/min_x); % linear distribution with minimum grid size
        vfine = ppval(cs,xfine);
        [vmax, pos] = max(vfine);
        vmax_list(j) = vmax;
        xmax_list(j) = xfine(pos); %1-xfine(pos);
        
        
        [vmin, pos] = min(vmid);
        disp(vmin);
        disp(xp(pos));   %    disp(1-xp(pos));
        [vmin, pos] = min(vfine);
        vmin_list(j) = vmin;
        xmin_list(j) = 1-xfine(pos);
        
        figure
        hold on
        %     if (j==1)
        %     plot(xBP,vBP,'bx');
        %     end
        %     plot(xBP,vGhia,'rx');
        plot(xp,v(:,Ny/2),color)
        xlabel('x')
        ylabel('v')
        box
        title('vertical velocity at horizontal centerline');
        
        figure
        hold on
        %     plot(xBP,whBP,'bx');
        plot(x(2:end-1),-omega(:,Ny/2),color)
        xlabel('x')
        ylabel('\omega')
        box
        title('vorticity at horizontal centerline');
        
        figure
        hold on
        %     if (j==1)
        %     plot(xBP,phBP,'bx');
        %     end
        px = (pres_(:,Ny/2)+pres_(:,Ny/2+1))/2;
        plot(xp,px,color)
        xlabel('x')
        ylabel('p')
        box
        title('pressure at horizontal centerline');
        
    else
        [vmax, pos] = max((v(:,ceil(Ny/2))+v(:,ceil(Ny/2)))/2);
        disp(vmax);
        disp(xp(pos)); %disp(1-xp(pos));
        [vmin, pos] = min((v(:,ceil(Ny/2))+v(:,ceil(Ny/2)))/2);
        disp(vmin);
        disp(xp(pos)); %disp(1-xp(pos));
    end
    
    
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
    
end

%% write output
if (write_files == 1 && save_results == 1)
    
    file_uy  = ['../../' path_results,'/u_y_t1_Re100.txt'];
    data = [yp u(Nx/2,:)'];
    % fid = fopen(file_uy,'wt');
    % fprintf(fid,'%.8f\n',data);
    % fclose(fid);
    save(file_uy, 'data', '-ascii', '-double', '-tabs')
    
    
    file_vx  = ['../../' path_results,'/v_x_t1_Re100.txt'];
    data = [xp v(:,Ny/2)];
    % fid = fopen(file_uy,'wt');
    % fprintf(fid,'%.8f\n',data);
    % fclose(fid);
    save(file_vx, 'data', '-ascii', '-double', '-tabs')
    
    file_py  = ['../../' path_results,'/p_y_t1_Re100.txt'];
    data = [yp py'];
    % fid = fopen(file_uy,'wt');
    % fprintf(fid,'%.8f\n',data);
    % fclose(fid);
    save(file_py, 'data', '-ascii', '-double', '-tabs')
    
    file_px  = ['../../' path_results,'/p_x_t1_Re100.txt'];
    data = [xp px];
    % fid = fopen(file_uy,'wt');
    % fprintf(fid,'%.8f\n',data);
    % fclose(fid);
    save(file_px, 'data', '-ascii', '-double', '-tabs')
    
    
end

