%% post-processing Lid-driven cavity

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));

%load BP data
% run('results/LDC/BP.m');

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);

if (floor(Nx/2)==Nx/2 && floor(Ny/2)==Ny/2)
    pres_ = pres-(pres(Nx/2+1,Ny/2+1)+pres(Nx/2,Ny/2))/2;
else
    pres_ = pres-pres(ceil(Nx/2),ceil(Ny/2));
end

%% create centerline plots

if (floor(Nx/2)==Nx/2) % Nx even
    disp('min u along vertical centerline');
    
    umid = u(Nx/2,:);
    [umax pos] = max(umid);
    disp(umax);
    disp(yp(pos));  
    
    %reconstruct data with spline for more accurate result:
    if (min(hy) == max(hy)) % uniform grid
        min_y = deltay/10;
    else
        min_y = min(hy);
    end
    cs = spline(yp,umid);
    yfine = linspace(yp(1),yp(end),L_y/min_y); % linear distribution with minimum grid size
    ufine = ppval(cs,yfine);
    [umax pos] = max(ufine);
    umax_list(j) = umax;
    ymax_list(j) = yfine(pos);
    
    %%
    figure(1)
    hold on
%     if (j==1)
%     plot(uBP,yBP,'bx');
%     end
%     plot(-uGhia,yBP,'rx');
    plot(u(Nx/2,:),yp,color) 
    xlabel('u')
    ylabel('y')
    box
    
    %%
%     figure(2)
%     hold on
% %     plot(wvBP,yBP,'bx');
%     plot(-omega(Nx/2,:),y(2:end-1),color) 
%     xlabel('w')
%     ylabel('y')
%     box
    
    figure(3)
    hold on
%     if (j==1)    
%     plot(pvBP,yBP,'bx');
%     end
    plot((pres_(Nx/2,:)+pres_(Nx/2+1,:))/2,yp,color) 
    xlabel('p')
    ylabel('y')
    box
else
    [umin pos] = min((u(ceil(Nx/2),:)+u(floor(Nx/2),:))/2);
    disp(umin);
    disp(yp(pos));
end

if (floor(Ny/2)==Ny/2) % Ny even
    disp('max and min v along vertical centerline');
    
    vmid  = v(:,Ny/2);      
    [vmax pos] = max(vmid);
    disp(vmax);
    disp(1-xp(pos));   
    
    %reconstruct data with spline for more accurate results:
    if (min(hx) == max(hx)) % uniform grid
        min_x = deltax/10;
    else
        min_x = min(hx);
    end    
    cs = spline(xp,vmid);
    xfine = linspace(xp(1),xp(end),L_x/min_x); % linear distribution with minimum grid size
    vfine = ppval(cs,xfine);
    [vmax pos] = max(vfine);
    vmax_list(j) = vmax;
    xmax_list(j) = 1-xfine(pos);
    
    
    [vmin pos] = min(vmid);
    disp(vmin);
    disp(1-xp(pos));       
    [vmin pos] = min(vfine);
    vmin_list(j) = vmin;
    xmin_list(j) = 1-xfine(pos);
   
    figure(4)
    hold on
%     if (j==1)
%     plot(xBP,vBP,'bx');
%     end
%     plot(xBP,vGhia,'rx');
    plot(xp,v(:,Ny/2),color) 
    xlabel('x')
    ylabel('v')
    box
%     title('vertical velocity at horizontal centerline');
    
%     figure(5)
%     hold on
% %     plot(xBP,whBP,'bx');
%     plot(x(2:end-1),-omega(:,Ny/2),color)
%     xlabel('x')
%     ylabel('\omega')
%     box
%     title('vorticity at horizontal centerline');

    figure(6)
    hold on
%     if (j==1)
%     plot(xBP,phBP,'bx');
%     end
    plot(xp,(pres_(:,Ny/2)+pres_(:,Ny/2+1))/2,color)
    xlabel('x')
    ylabel('p')
    box
%     title('pressure at horizontal centerline');

else
    [vmax pos] = max((v(:,ceil(Ny/2))+v(:,ceil(Ny/2)))/2);
    disp(vmax);
    disp(1-xp(pos));
    [vmin pos] = min((v(:,ceil(Ny/2))+v(:,ceil(Ny/2)))/2);
    disp(vmin);
    disp(1-xp(pos));
end


%% convergence graphs

if (j==length(mesh_list) && length(mesh_list)>1)    
    
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