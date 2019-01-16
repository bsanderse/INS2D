
%% pressure volumes

% number of pressure points
Npx         = Nx;
Npy         = Ny;
Np          = Npx*Npy;


%% u-volumes

% x(1)   x(2)   x(3) ....      x(Nx)   x(Nx+1)
% |      |      |              |       |
% |      |      |              |       |
% Dirichlet BC: 
% uLe    u(1)   u(2) ....      u(Nx-1) uRi
% periodic BC: 
% u(1)   u(2)   u(3) ....      u(Nx)   u(1)
% pressure BC:
% u(1)   u(2)   u(3) ....      u(Nx)   u(Nx+1) 


% x-dir
Nux_b       = 2;               % boundary points
Nux_in      = Nx+1;            % inner points
if (strcmp(BC.u.left,'dir') || strcmp(BC.u.left,'sym'))
    Nux_in  = Nux_in-1;
end
if (strcmp(BC.u.right,'dir') || strcmp(BC.u.right,'sym'))
    Nux_in  = Nux_in-1;
end
if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
    Nux_in  = Nux_in-1;
end
Nux_t       = Nux_in + Nux_b;  % total number 

% y-dir
Nuy_b       = 2;               % boundary points
Nuy_in      = Ny;              % inner points
Nuy_t       = Nuy_in + Nuy_b;  % total number

% total number 
Nu          = Nux_in*Nuy_in;  


%% v-volumes

% x-dir
Nvx_b       = 2;               % boundary points
Nvx_in      = Nx;              % inner points
Nvx_t       = Nvx_in + Nvx_b;  % total number

% y-dir
Nvy_b       = 2;               % boundary points
Nvy_in      = Ny+1;            % inner points
if (strcmp(BC.v.low,'dir') || strcmp(BC.v.low,'sym'))
    Nvy_in  = Nvy_in-1;
end
if (strcmp(BC.v.up,'dir') || strcmp(BC.v.up,'sym'))
    Nvy_in  = Nvy_in-1;
end
if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
    Nvy_in  = Nvy_in-1;
end
Nvy_t       = Nvy_in + Nvy_b;  % total number 

% total number 
Nv          = Nvx_in*Nvy_in;  


%% extra variables 
N1 = (Nux_in+1)*Nuy_in; %size(Iu_ux,1);
N2 = Nux_in*(Nuy_in+1); %size(Iv_uy,1);
N3 = (Nvx_in+1)*Nvy_in; % size(Iu_vx,1); 
N4 =  Nvx_in*(Nvy_in+1); % size(Iv_vy,1);


%% for a grid with three times larger volumes:
if (order4==1)
    hx3          = zeros(Nx,1);
    hx3(2:end-1) = hx(1:end-2)+hx(2:end-1)+hx(3:end);
    if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
        hx3(1)   = hx(end)+hx(1)+hx(2);
        hx3(end) = hx(end-1)+hx(end)+hx(1);
    else
        hx3(1)   = 2*hx(1)+hx(2);
        hx3(end) = hx(end-1)+2*hx(end);
    end

    hy3          = zeros(Ny,1);
    hy3(2:end-1) = hy(1:end-2)+hy(2:end-1)+hy(3:end);
    if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
        hy3(1)   = hy(end)+hy(1)+hy(2);
        hy3(end) = hy(end-1)+hy(end)+hy(1);
    else
        hy3(1)   = 2*hy(1)+hy(2);
        hy3(end) = hy(end-1)+2*hy(end);
    end    

    % distance between pressure points
    gx3          = zeros(Nx+1,1);
    gx3(3:Nx-1)  = gx(2:end-3)+gx(3:end-2)+gx(4:end-1);
    if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
        gx3(1)   = gx(end-1)+gx(end)+gx(1)+gx(2);
        gx3(2)   = gx(end)+gx(1)+gx(2)+gx(3);
        gx3(end-1) = gx(end-2)+gx(end-1)+gx(end)+gx(1);
        gx3(end) = gx(end-1)+gx(end)+gx(1)+gx(2);
    else
        gx3(1)   = 2*gx(1)+2*gx(2);
        gx3(2)   = 2*gx(1)+gx(2)+gx(3);
        gx3(end-1) = 2*gx(end)+gx(end-1)+gx(end-2); 
        gx3(end) = 2*gx(end)+2*gx(end-1);
    end

    % distance between pressure points
    gy3          = zeros(Ny+1,1);
    gy3(3:Ny-1)  = gy(2:end-3)+gy(3:end-2)+gx(4:end-1);
    if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
        gy3(1)   = gy(end-1)+gy(end)+gy(1)+gy(2);
        gy3(2)   = gy(end)+gy(1)+gy(2)+gy(3);
        gy3(end-1) = gy(end-2)+gy(end-1)+gy(end)+gy(1);
        gy3(end) = gy(end-1)+gy(end)+gy(1)+gy(2);
    else
        gy3(1)   = 2*gy(1)+2*gy(2);
        gy3(2)   = 2*gy(1)+gy(2)+gy(3);
        gy3(end-1) = 2*gy(end)+gy(end-1)+gy(end-2); 
        gy3(end) = 2*gy(end)+2*gy(end-1);
    end    

end

%% adapt mesh metrics depending on number of volumes

%% x-direction

% gxi: integration and gxd: differentiation
gxd         = gx;
gxd(1)      = hx(1);
gxd(end)    = hx(end);

% hxi: integration and hxd: differentiation
% map to find suitable size
hxi         = hx;
hxd         = [hx(1); hx; hx(end)];

% restrict Nx+2 to Nux_in+1 points
if (strcmp(BC.u.left,'dir') && strcmp(BC.u.right,'dir'))
    diagpos = 1;
end
if (strcmp(BC.u.left,'dir') && strcmp(BC.u.right,'pres'))
    diagpos = 1;
end
if (strcmp(BC.u.left,'pres') && strcmp(BC.u.right,'dir'))
    diagpos = 0;
end
if (strcmp(BC.u.left,'pres') && strcmp(BC.u.right,'pres'))
    diagpos = 0; 
end
if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
    diagpos = 0;
end
Bmap  = spdiags(ones(Nux_in+1,1),diagpos,Nux_in+1,Nx+2);
Bmap2 = spdiags(ones(Nux_in,1),diagpos,Nux_in,Nx+1);
hxd   = Bmap*hxd;
gxi   = Bmap2*gx;

if (order4==1)
    hxi3  = hx3;
    
%     if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))  
    hxd3  = [hx3(end-1); hx3(end); hx3; hx3(1); hx3(2)];
    gxd3  = [2*gx(2)+2*gx(1); 2*gx(1)+2*gx(2); gx3(2:end-1); 2*gx(end)+2*gx(end-1); 2*gx(end-1)+2*gx(end)];
    % map
    Bmap3 = spdiags(ones(Nux_in+3,1),diagpos,Nux_in+3,Nx+4);   
    gxi3  = Bmap2*gx3;
    hxd3  = Bmap3*hxd3;
    
    hxd13  = [hxd(1);hxd;hxd(end)];
    gyd13  = [gyd(1);gyd;gyd(end)];
    gxd13  = [gxd(1);gxd;gxd(end)];
    hyd13  = [hyd(1);hyd;hyd(end)];
    
    if (strcmp(BC.u.left,'pres') || strcmp(BC.u.right,'pres'))
        error('not implemented');
    end
    
    
end


% coordinates of velocity points with BC taken into account:
xin   = Bmap2*x;  
% y coordinate is yp

% adapt some metrics in case of periodic BC
if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
    gxd(1)  = (hx(1)+hx(end))/2;
    gxd(end)= gxd(1);    
    gxi(1)  = gxd(1);
    hxd(1)  = hx(end);
end

% matrix to map from Nvx_t-1 to Nux_in points 
% (used in interpolation, convection_diffusion, viscosity)
Bvux  = spdiags(ones(Nvx_t-1,1),diagpos,Nux_in,Nvx_t-1);
% if (order4==1)
%     Bvux3 = spdiags(ones(Nvx_t-1,1),diagpost,Nux_in+2,Nvx_t+1);
% end
% map from Npx+2 points to Nux_t-1 points (ux faces)
Bkux  = Bmap;

%% y-direction

% gyi: integration and gyd: differentiation
gyd         = gy;
gyd(1)      = hy(1);
gyd(end)    = hy(end);

% hyi: integration and hyd: differentiation
% map to find suitable size
hyi         = hy;
hyd         = [hy(1); hy; hy(end)];

if (order4==1)
    hyi3        = hy3;
    hyd3        = [hy3(1); hy3(1); hy3; hy3(end); hy3(end)];
    gyd3        = [2*gy(2)+2*gy(1); 2*gy(1)+2*gy(2); gy3(2:end-1); 2*gy(end)+2*gy(end-1); 2*gy(end-1)+2*gy(end)];
end

% restrict Ny+2 to Nvy_in+1 points
if (strcmp(BC.v.low,'dir') && strcmp(BC.v.up,'dir'))
    diagpos = 1;
end
if (strcmp(BC.v.low,'dir') && strcmp(BC.v.up,'pres'))
    diagpos = 1;
end
if (strcmp(BC.v.low,'pres') && strcmp(BC.v.up,'dir'))
    diagpos = 0;
end
if (strcmp(BC.v.low,'pres') && strcmp(BC.v.up,'pres'))
    diagpos = 0; 
end
if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
    diagpos = 0;
end
Bmap  = spdiags(ones(Nvy_in+1,1),diagpos,Nvy_in+1,Ny+2);
Bmap2 = spdiags(ones(Nvy_in,1),diagpos,Nvy_in,Ny+1);

hyd   = Bmap*hyd;
gyi   = Bmap2*gy;

if (order4==1)
    Bmap3 = spdiags(ones(Nvy_in+3,1),diagpos,Nvy_in+3,Ny+4);    
    gyi3  = Bmap2*gy3;
    hyd3  = Bmap3*hyd3;
end


% coordinates of velocity points with BC taken into account
yin   = Bmap2*y;  
% x coordinate is xp

% adapt some metrics in case of periodic BC
if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
    gyd(1)  = (hy(1)+hy(end))/2;
    gyd(end)= gyd(1);    
    gyi(1)  = gyd(1);
    hyd(1)  = hy(end);
end

% matrix to map from Nuy_t-1 to Nvy_in points 
% (used in interpolation, convection_diffusion)
Buvy   = spdiags(ones(Nuy_t-1,1),diagpos,Nvy_in,Nuy_t-1);
% map from Npy+2 points to Nvy_t-1 points (vy faces)
Bkvy  = Bmap;

% volume (area) of pressure control volumes
Omp    = kron(hyi,hxi);
Omp_inv= 1./Omp;
% volume (area) of u control volumes
Omu    = kron(hyi,gxi);
Omu_inv= 1./Omu;
% volume of ux volumes
Omux   = kron(hyi,hxd);
% volume of uy volumes
Omuy   = kron(gyd,gxi);
% volume (area) of v control volumes
Omv    = kron(gyi,hxi);
Omv_inv= 1./Omv;
% volume of vx volumes
Omvx   = kron(gyi,gxd);
% volume of vy volumes
Omvy   = kron(hyd,hxi);
% volume (area) of vorticity control volumes
Omvort = kron(gyi,gxi);
Omvort_inv = 1./Omvort;

Om     = [Omu; Omv];
Om_inv = [Omu_inv; Omv_inv];


if (order4==1)
    
    % differencing for second order operators on the fourth order mesh
    Omux1  = kron(hyi,hxd13);
    Omuy1  = kron(gyd13,gxi);
    Omvx1  = kron(gyi,gxd13);
    Omvy1  = kron(hyd13,hxi);
    
    % volume (area) of pressure control volumes
    Omp3  = kron(hyi3,hxi3); 
    % volume (area) of u-vel control volumes
    Omu3  = kron(hyi3,gxi3);
    % volume (area) of v-vel control volumes
    Omv3  = kron(gyi3,hxi3);
    % volume (area) of dudx control volumes
    Omux3 = kron(hyi3,hxd3);
    % volume (area) of dudy control volumes
    Omuy3 = kron(gyd3,gxi3);
    % volume (area) of dvdx control volumes
    Omvx3 = kron(gyi3,gxd3);
    % volume (area) of dvdy control volumes
    Omvy3 = kron(hyd3,hxi3);
    
    Omu1  = Omu;
    Omv1  = Omv;
    
    Omu   = alfa*Omu1 - Omu3;
    Omv   = alfa*Omv1 - Omv3;
    Omu_inv = 1./Omu;
    Omv_inv = 1./Omv;
    Om    = [Omu;Omv];
    Om_inv= [Omu_inv; Omv_inv];
    
    Omux  = alfa*Omux1 - Omux3;
    Omuy  = alfa*Omuy1 - Omuy3;
    Omvx  = alfa*Omvx1 - Omvx3;
    Omvy  = alfa*Omvy1 - Omvy3;
    
end


% metrics that can be useful for initialization:
xu = kron(ones(1,Nuy_in),xin);
yu = kron(yp,ones(Nux_in,1));
xu = reshape(xu,Nux_in,Nuy_in);
yu = reshape(yu,Nux_in,Nuy_in);

xv = kron(ones(1,Nvy_in),xp);
yv = kron(yin,ones(Nvx_in,1));
xv = reshape(xv,Nvx_in,Nvy_in);
yv = reshape(yv,Nvx_in,Nvy_in);

xpp = kron(ones(Ny,1),xp);
ypp = kron(yp,ones(Nx,1));
xpp = reshape(xpp,Nx,Ny);
ypp = reshape(ypp,Nx,Ny);
  




% plot the grid: velocity points and pressure points
if plotgrid==1
    plot_staggered(x,y);    
end