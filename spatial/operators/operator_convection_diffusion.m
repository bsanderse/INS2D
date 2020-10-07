function options = operator_convection_diffusion(options)
% construct convection and diffusion operators

% boundary conditions
BC = options.BC;

Nx = options.grid.Nx;
Ny = options.grid.Ny;

% number of interior points and boundary points
Nu     = options.grid.Nu;
Nux_in = options.grid.Nux_in;
Nux_b  = options.grid.Nux_b;
Nux_t  = options.grid.Nux_t;
Nuy_in = options.grid.Nuy_in;
Nuy_b  = options.grid.Nuy_b;
Nuy_t  = options.grid.Nuy_t;
Nv     = options.grid.Nv;
Nvx_in = options.grid.Nvx_in;
Nvx_b  = options.grid.Nvx_b;
Nvx_t  = options.grid.Nvx_t;
Nvy_in = options.grid.Nvy_in;
Nvy_b  = options.grid.Nvy_b;
Nvy_t  = options.grid.Nvy_t;

hx  = options.grid.hx;
hy  = options.grid.hy;
hxi = options.grid.hxi;
hyi = options.grid.hyi;
hxd = options.grid.hxd;
hyd = options.grid.hyd;
gxi = options.grid.gxi;
gyi = options.grid.gyi;
gxd = options.grid.gxd;
gyd = options.grid.gyd;

Buvy = options.grid.Buvy;
Bvux = options.grid.Bvux;

order4 = options.discretization.order4;

if (order4 == 1)
    alfa   = options.discretization.alfa;
    hxi3   = options.grid.hxi3;
    hyi3   = options.grid.hyi3;
    gxi3   = options.grid.gxi3;
    gyi3   = options.grid.gyi3;
    hxd13  = options.grid.hxd13;
    hxd3   = options.grid.hxd3;
    hyd13  = options.grid.hyd13;
    hyd3   = options.grid.hyd3;
    gxd13  = options.grid.gxd13;
    gxd3   = options.grid.gxd3;
    gyd13  = options.grid.gyd13;
    gyd3   = options.grid.gyd3;
    Omux   = options.grid.Omux;
    Omuy   = options.grid.Omuy;
    Omvx   = options.grid.Omvx;
    Omvy   = options.grid.Omvy;
    Omux1  = options.grid.Omux1;
    Omux3  = options.grid.Omux3;
    Omuy1  = options.grid.Omuy1;
    Omuy3  = options.grid.Omuy3;
    Omvx1  = options.grid.Omvx1;
    Omvx3  = options.grid.Omvx3;
    Omvy1  = options.grid.Omvy1;
    Omvy3  = options.grid.Omvy3;
end

visc = options.case.visc;
Re   = options.fluid.Re;


%% Convection (differencing) operator Cu

% calculates difference from pressure points to velocity points
diag1       = ones(Nux_t,1);
D1D         = spdiags([-diag1 diag1],[0 1],Nux_t-2,Nux_t-1);
Cux         = kron(speye(Nuy_in),D1D);
if (order4==0)
    Dux         = kron(spdiags(hyi,0,Ny,Ny),D1D);
end
clear diag1 D1D

% calculates difference from corner points to velocity points
diag1       = ones(Nuy_t,1);
D1D         = spdiags([-diag1 diag1],[0 1],Nuy_t-2,Nuy_t-1);
Cuy         = kron(D1D,speye(Nux_in));
if (order4==0)
    Duy         = kron(D1D,spdiags(gxi,0,Nux_in,Nux_in));
end
clear diag1 D1D

% Cu          = [Cux Cuy];
% Du          = [Dux Duy];


%% Convection (differencing) operator Cv

% calculates difference from pressure points to velocity points
diag1       = ones(Nvx_t,1);
D1D         = spdiags([-diag1 diag1],[0 1],Nvx_t-2,Nvx_t-1);
Cvx         = kron(speye(Nvy_in),D1D);
if (order4==0)
    Dvx         = kron(spdiags(gyi,0,Nvy_in,Nvy_in),D1D);
end
clear diag1 D1D

% calculates difference from corner points to velocity points
diag1       = ones(Nvy_t,1);
D1D         = spdiags([-diag1 diag1],[0 1],Nvy_t-2,Nvy_t-1);
Cvy         = kron(D1D,speye(Nvx_in));
if (order4==0)
    Dvy         = kron(D1D,spdiags(hxi,0,Nx,Nx));
end
clear diag1 D1D

% Cv          = [Cvx Cvy];
% Dv          = [Dvx Dvy];

if (order4==0)
    %% Diffusion operator (stress tensor), u-component
    % similar to averaging, but with mesh sizes
    
    %% Su_ux: evaluate ux
    diag1        = 1./hxd;
    S1D          = spdiags([-diag1 diag1],[0 1],Nux_t-1,Nux_t);
    
    % boundary conditions
    Su_ux_BC     = BC_general(Nux_t,Nux_in,Nux_b, ...
        BC.u.left,BC.u.right,hx(1),hx(end));
    
    % extend to 2D
    Su_ux        = kron(speye(Ny),S1D*Su_ux_BC.B1D);
    Su_ux_BC.Bbc = kron(speye(Ny),S1D*Su_ux_BC.Btemp);
    
    clear diag1 S1D
    
    %% Su_uy: evaluate uy
    diag1        = 1./gyd;
    S1D          = spdiags([-diag1 diag1],[0 1],Nuy_t-1,Nuy_t);
    
    % boundary conditions
    % Su_uy_BC     = BC_general_stag(Nuy_t,Nuy_in,Nuy_b, ...
    %                                            BC.u.low,BC.u.up,hy(1),hy(end));
    Su_uy_BC     = BC_diff_stag(Nuy_t,Nuy_in,Nuy_b, ...
        BC.u.low,BC.u.up,hy(1),hy(end));
    
    
    % extend to 2D
    Su_uy        = kron(S1D*Su_uy_BC.B1D,speye(Nux_in));
    Su_uy_BC.Bbc = kron(S1D*Su_uy_BC.Btemp,speye(Nux_in));
    
    clear diag1 S1D
    
    %% Sv_uy: evaluate vx at uy;
    % same as Iv_uy except for mesh sizes and -diag diag
    
    diag1           = 1./gxd;
    S1D             = spdiags([-diag1 diag1],[0 1],Nvx_t-1,Nvx_t);
    % the restriction is essentially 1D so it can be directly applied to I1D
    S1D             = Bvux*S1D;
    S2D             = kron(speye(Nuy_t-1),S1D);
    
    
    % boundary conditions low/up
    Nb              = Nuy_in+1-Nvy_in;
    Sv_uy_BC_lu     = BC_general(Nuy_in+1,Nvy_in,Nb,...
        BC.v.low,BC.v.up,hy(1),hy(end));
    Sv_uy_BC_lu.B2D = kron(Sv_uy_BC_lu.B1D,speye(Nvx_in));
    Sv_uy_BC_lu.Bbc = kron(Sv_uy_BC_lu.Btemp,speye(Nvx_in));
    
    
    % boundary conditions left/right
    Sv_uy_BC_lr     = BC_general_stag(Nvx_t,Nvx_in,Nvx_b,...
        BC.v.left,BC.v.right,hx(1),hx(end));
    % take I2D into left/right operators for convenience
    Sv_uy_BC_lr.B2D = S2D*kron(speye(Nuy_t-1),Sv_uy_BC_lr.B1D);
    Sv_uy_BC_lr.Bbc = S2D*kron(speye(Nuy_t-1),Sv_uy_BC_lr.Btemp);
    
    
    % resulting operator:
    Sv_uy           = Sv_uy_BC_lr.B2D*Sv_uy_BC_lu.B2D;
    
    clear diag1 S1D
    
    %% Diffusion operator (stress tensor), v-component
    % similar to averaging!
    
    %% Su_vx: evaluate uy at vx
    % same as Iu_vx except for mesh sizes and -diag diag
    
    diag1           = 1./gyd;
    S1D             = spdiags([-diag1 diag1],[0 1],Nuy_t-1,Nuy_t);
    S1D             = Buvy*S1D;
    S2D             = kron(S1D,speye(Nvx_t-1));
    
    % boundary conditions low/up
    Su_vx_BC_lu     = BC_general_stag(Nuy_t,Nuy_in,Nuy_b,...
        BC.u.low,BC.u.up,hy(1),hy(end));
    Su_vx_BC_lu.B2D = S2D*kron(Su_vx_BC_lu.B1D,speye(Nvx_t-1));
    Su_vx_BC_lu.Bbc = S2D*kron(Su_vx_BC_lu.Btemp,speye(Nvx_t-1));
    
    % boundary conditions left/right
    Nb              = Nvx_in+1-Nux_in;
    Su_vx_BC_lr     = BC_general(Nvx_in+1,Nux_in,Nb,...
        BC.u.left,BC.u.right,hx(1),hx(end));
    
    Su_vx_BC_lr.B2D = kron(speye(Nuy_in),Su_vx_BC_lr.B1D);
    Su_vx_BC_lr.Bbc = kron(speye(Nuy_in),Su_vx_BC_lr.Btemp);
    
    % resulting operator:
    Su_vx           = Su_vx_BC_lu.B2D*Su_vx_BC_lr.B2D;
    
    clear diag1 S1D
    
    %% Sv_vx: evaluate vx
    diag1        = 1./gxd;
    S1D          = spdiags([-diag1 diag1],[0 1],Nvx_t-1,Nvx_t);
    
    % boundary conditions
    % Sv_vx_BC     = BC_general_stag(Nvx_t,Nvx_in,Nvx_b,...
    %                                            BC.v.left,BC.v.right,hx(1),hx(end));
    Sv_vx_BC     = BC_diff_stag(Nvx_t,Nvx_in,Nvx_b,...
        BC.v.left,BC.v.right,hx(1),hx(end));
    
    % extend to 2D
    Sv_vx        = kron(speye(Nvy_in),S1D*Sv_vx_BC.B1D);
    Sv_vx_BC.Bbc = kron(speye(Nvy_in),S1D*Sv_vx_BC.Btemp);
    
    clear diag1 S1D
    
    %% Sv_vy: evaluate vy
    diag1        = 1./hyd;
    S1D          = spdiags([-diag1 diag1],[0 1],Nvy_t-1,Nvy_t);
    
    % boundary conditions
    Sv_vy_BC     = BC_general(Nvy_t,Nvy_in,Nvy_b, ...
        BC.v.low,BC.v.up,hy(1),hy(end));
    
    % extend to 2D
    Sv_vy        = kron(S1D*Sv_vy_BC.B1D,speye(Nx));
    Sv_vy_BC.Bbc = kron(S1D*Sv_vy_BC.Btemp,speye(Nx));
    
    clear diag1 S1D
    
end

%% fourth order operators
if (order4==1)
    
    %% Convection (differencing) operator Cu
    
    % calculates difference from pressure points to velocity points
    diag1       = ones(Nux_t,1);
    D1D         = spdiags([-diag1 diag1],[1 2],Nux_t-2,Nux_t+1);
    Dux         = kron(spdiags(hyi,0,Ny,Ny),D1D);
    % the 'second order' Cux is unchanged
    % the 'second order' Dux changes, because we also use the 'second
    % order' flux at 'fourth order' ghost points (Dux should have the same
    % size as Dux3)
    
    clear diag1 D1D
    
    % calculates difference from pressure points to velocity points
    diag1       = ones(Nux_t,1);
    D1D3        = spdiags([-diag1 diag1],[0 3],Nux_t-2,Nux_t+1);
    Cux3        = kron(speye(Ny),D1D3);
    Dux3        = kron(spdiags(hyi3,0,Ny,Ny),D1D3);
    
    clear diag1 D1D3
    
    % calculates difference from corner points to velocity points
    diag1       = ones(Nuy_t,1);
    D1D         = spdiags([-diag1 diag1],[1 2],Nuy_t-2,Nuy_t+1);
    Duy         = kron(D1D,spdiags(gxi,0,Nux_in,Nux_in));
    
    clear diag1 D1D
    
    % calculates difference from corner points to velocity points
    diag1       = ones(Nuy_t,1);
    D1D3        = spdiags([-diag1 diag1],[0 3],Nuy_t-2,Nuy_t+1);
    % uncomment for new BC (functions/new)
    if ( strcmp(BC.u.low,'dir') )
        D1D3(1,1) = 1; D1D3(1,2) = -2;
    end
    if ( strcmp(BC.u.up,'dir') )
        D1D3(end,end-1) = 2; D1D3(end,end) = -1;
    end
    Cuy3        = kron(D1D3,speye(Nux_in));
    Duy3        = kron(D1D3,spdiags(gxi3,0,Nux_in,Nux_in));
    
    clear diag1 D1D3
    
    
    %% Convection (differencing) operator Cv
    
    % calculates difference from pressure points to velocity points
    diag1       = ones(Nvx_t,1);
    D1D         = spdiags([-diag1 diag1],[1 2],Nvx_t-2,Nvx_t+1);
    Dvx         = kron(spdiags(gyi,0,Nvy_in,Nvy_in),D1D);
    
    clear diag1 D1D
    
    % calculates difference from pressure points to velocity points
    diag1       = ones(Nvx_t,1);
    D1D3        = spdiags([-diag1 diag1],[0 3],Nvx_t-2,Nvx_t+1);
    % uncomment for new BC (functions/new)
    if ( strcmp(BC.v.left,'dir') )
        D1D3(1,1) = 1; D1D3(1,2) = -2;
    end
    if ( strcmp(BC.v.right,'dir') )
        D1D3(end,end-1) = 2; D1D3(end,end) = -1;
    end
    Cvx3        = kron(speye(Nvy_in),D1D3);
    Dvx3        = kron(spdiags(gyi3,0,Nvy_in,Nvy_in),D1D3);
    
    clear diag1 D1D3
    
    % calculates difference from corner points to velocity points
    diag1       = ones(Nvy_t,1);
    D1D         = spdiags([-diag1 diag1],[1 2],Nvy_t-2,Nvy_t+1);
    Dvy         = kron(D1D,spdiags(hxi,0,Nx,Nx));
    
    clear diag1 D1D
    
    % calculates difference from corner points to velocity points
    diag1       = ones(Nvy_t,1);
    D1D3        = spdiags([-diag1 diag1],[0 3],Nvy_t-2,Nvy_t+1);
    Cvy3        = kron(D1D3,speye(Nvx_in));
    Dvy3        = kron(D1D3,spdiags(hxi3,0,Nx,Nx));
    
    clear diag1 D1D3
    
    
    %% Su_ux: evaluate ux
    diag1        = 1./hxd13;
    S1D          = spdiags([-diag1 diag1],[1 2],Nux_in+3,Nux_t+4);
    
    % boundary conditions
    Su_ux_BC     = BC_diff3(Nux_t+4,Nux_in,Nux_t+4-Nux_in, ...
        BC.u.left,BC.u.right,hx(1),hx(end));
    
    % extend to 2D
    Su_ux        = spdiags(Omux1,0,length(Omux1),length(Omux1))*kron(speye(Ny),S1D*Su_ux_BC.B1D);
    Su_ux_BC.Bbc = spdiags(Omux1,0,length(Omux1),length(Omux1))*kron(speye(Ny),S1D*Su_ux_BC.Btemp);
    
    clear diag1 S1D
    
    diag1        = 1./hxd3;
    S1D3         = spdiags([-diag1 diag1],[0 3],Nux_in+3,Nux_t+4);
    
    % boundary conditions
    Su_ux_BC3    = BC_diff3(Nux_t+4,Nux_in,Nux_t+4-Nux_in, ...
        BC.u.left,BC.u.right,hx(1),hx(end));
    % extend to 2D
    Su_ux3        = spdiags(Omux3,0,length(Omux3),length(Omux3))*kron(speye(Nuy_in),S1D3*Su_ux_BC3.B1D);
    Su_ux_BC3.Bbc = spdiags(Omux3,0,length(Omux3),length(Omux3))*kron(speye(Nuy_in),S1D3*Su_ux_BC3.Btemp);
    
    clear diag1 S1D3
    
    
    %% Su_uy: evaluate uy
    diag1        = 1./gyd13;
    S1D          = spdiags([-diag1 diag1],[1 2],Nuy_in+3,Nuy_t+4);
    % boundary conditions
    Su_uy_BC     = BC_diff_stag3(Nuy_t+4,Nuy_in,Nuy_t+4-Nuy_in, ...
        BC.u.low,BC.u.up,hy(1),hy(end));
    % extend to 2D
    Su_uy        = spdiags(Omuy1,0,length(Omuy1),length(Omuy1))*kron(S1D*Su_uy_BC.B1D,speye(Nux_in));
    Su_uy_BC.Bbc = spdiags(Omuy1,0,length(Omuy1),length(Omuy1))*kron(S1D*Su_uy_BC.Btemp,speye(Nux_in));
    clear diag1 S1D
    
    
    diag1        = 1./gyd3;
    S1D3         = spdiags([-diag1 diag1],[0 3],Nuy_in+3,Nuy_t+4);
    % boundary conditions
    Su_uy_BC3    = BC_diff_stag3(Nuy_t+4,Nuy_in,Nuy_t+4-Nuy_in, ...
        BC.u.low,BC.u.up,hy(1),hy(end));
    % extend to 2D
    Su_uy3        = spdiags(Omuy3,0,length(Omuy3),length(Omuy3))*kron(S1D3*Su_uy_BC3.B1D,speye(Nux_in));
    Su_uy_BC3.Bbc = spdiags(Omuy3,0,length(Omuy3),length(Omuy3))*kron(S1D3*Su_uy_BC3.Btemp,speye(Nux_in));
    
    clear diag1 S1D3
    
    
    %% Sv_vx: evaluate vx
    diag1        = 1./gxd13;
    S1D          = spdiags([-diag1 diag1],[1 2],Nvx_in+3,Nvx_t+4);
    
    % boundary conditions
    Sv_vx_BC     = BC_diff_stag3(Nvx_t+4,Nvx_in,Nvx_t+4-Nvx_in,...
        BC.v.left,BC.v.right,hx(1),hx(end));
    
    % extend to 2D
    Sv_vx        = spdiags(Omvx1,0,length(Omvx1),length(Omvx1))*kron(speye(Nvy_in),S1D*Sv_vx_BC.B1D);
    Sv_vx_BC.Bbc = spdiags(Omvx1,0,length(Omvx1),length(Omvx1))*kron(speye(Nvy_in),S1D*Sv_vx_BC.Btemp);
    
    clear diag1 S1D
    
    diag1        = 1./gxd3;
    S1D3         = spdiags([-diag1 diag1],[0 3],Nvx_in+3,Nvx_t+4);
    
    % boundary conditions
    Sv_vx_BC3    = BC_diff_stag3(Nvx_t+4,Nvx_in,Nvx_t+4-Nvx_in, ...
        BC.v.left,BC.v.right,hx(1),hx(end));
    % extend to 2D
    Sv_vx3        = spdiags(Omvx3,0,length(Omvx3),length(Omvx3))*kron(speye(Nvy_in),S1D3*Sv_vx_BC3.B1D);
    Sv_vx_BC3.Bbc = spdiags(Omvx3,0,length(Omvx3),length(Omvx3))*kron(speye(Nvy_in),S1D3*Sv_vx_BC3.Btemp);
    
    clear diag1 S1D3
    
    
    %% Sv_vy: evaluate vy
    diag1        = 1./hyd13;
    S1D          = spdiags([-diag1 diag1],[1 2],Nvy_in+3,Nvy_t+4);
    
    % boundary conditions
    Sv_vy_BC     = BC_diff3(Nvy_t+4,Nvy_in,Nvy_t+4-Nvy_in, ...
        BC.v.low,BC.v.up,hy(1),hy(end));
    
    % extend to 2D
    Sv_vy        = spdiags(Omvy1,0,length(Omvy1),length(Omvy1))*kron(S1D*Sv_vy_BC.B1D,speye(Nvx_in));
    Sv_vy_BC.Bbc = spdiags(Omvy1,0,length(Omvy1),length(Omvy1))*kron(S1D*Sv_vy_BC.Btemp,speye(Nvx_in));
    
    clear diag1 S1D
    
    diag1        = 1./hyd3;
    S1D3         = spdiags([-diag1 diag1],[0 3],Nvy_in+3,Nvy_t+4);
    
    % boundary conditions
    % Su_uy_BC     = BC_general_stag(Nuy_t,Nuy_in,Nuy_b, ...
    %                                            BC.u.low,BC.u.up,hy(1),hy(end));
    Sv_vy_BC3    = BC_diff3(Nvy_t+4,Nvy_in,Nvy_t+4-Nvy_in, ...
        BC.v.low,BC.v.up,hy(1),hy(end));
    % extend to 2D
    Sv_vy3        = spdiags(Omvy3,0,length(Omvy3),length(Omvy3))*kron(S1D3*Sv_vy_BC3.B1D,speye(Nvx_in));
    Sv_vy_BC3.Bbc = spdiags(Omvy3,0,length(Omvy3),length(Omvy3))*kron(S1D3*Sv_vy_BC3.Btemp,speye(Nvx_in));
    
    clear diag1 S1D3
    
end


%% assemble operators

switch visc
    
    case 'laminar'
        
        
        if (order4==0)
            Diffu  = Dux*( (1/Re)* Su_ux) + Duy*( (1/Re)* Su_uy);
            Diffv  = Dvx*( (1/Re)* Sv_vx) + Dvy*( (1/Re)* Sv_vy);
        elseif (order4==1)
            Diffux_div = (alfa*Dux - Dux3)*spdiags(1./Omux,0,length(Omux),length(Omux));
            Diffuy_div = (alfa*Duy - Duy3)*spdiags(1./Omuy,0,length(Omuy),length(Omuy));
            Diffvx_div = (alfa*Dvx - Dvx3)*spdiags(1./Omvx,0,length(Omvx),length(Omvx));
            Diffvy_div = (alfa*Dvy - Dvy3)*spdiags(1./Omvy,0,length(Omvy),length(Omvy));
            Diffu  = (1/Re)*Diffux_div*(alfa*Su_ux-Su_ux3) + ...
                     (1/Re)*Diffuy_div*(alfa*Su_uy-Su_uy3);
            Diffv  = (1/Re)*Diffvx_div*(alfa*Sv_vx-Sv_vx3) + ...
                     (1/Re)*Diffvy_div*(alfa*Sv_vy-Sv_vy3);
            
        end
        
        
    case {'keps','LES','qr','ML'}
        % only implemented for 2nd order
        
        % the terms below are an example of how the laminar case is
        % evaluated with the full stress tensor
        % these are not used in practical computations, as in the turbulent
        % case one needs to add nu_T, making the effective operator
        % solution-dependent, so that it cannot be computed beforehand
        
        % see diffusion.m for actual use
        
        % diffusion u-momentum
%         Diffu_u = Dux*( (1/Re) * 2*Su_ux) + Duy*( (1/Re) * Su_uy);
%         Diffu_v = Duy*( (1/Re) * Sv_uy);
        % diffusion v-momentum
%         Diffv_u = Dvx*( (1/Re) * Su_vx);
%         Diffv_v = Dvx*( (1/Re) * Sv_vx) + Dvy*( (1/Re) * 2*Sv_vy);
        
    otherwise
        error('wrong visc parameter');
end


options.discretization.Cux   = Cux;
options.discretization.Cuy   = Cuy;
options.discretization.Cvx   = Cvx;
options.discretization.Cvy   = Cvy;
options.discretization.Su_ux = Su_ux;
options.discretization.Su_uy = Su_uy;
options.discretization.Sv_vx = Sv_vx;
options.discretization.Sv_vy = Sv_vy;
options.discretization.Su_ux_BC = Su_ux_BC;
options.discretization.Su_uy_BC = Su_uy_BC;
options.discretization.Sv_vx_BC = Sv_vx_BC;
options.discretization.Sv_vy_BC = Sv_vy_BC;
options.discretization.Dux   = Dux;
options.discretization.Duy   = Duy;
options.discretization.Dvx   = Dvx;
options.discretization.Dvy   = Dvy;


switch visc

    case 'laminar'
        options.discretization.Diffu = Diffu;
        options.discretization.Diffv = Diffv;

    case {'keps','LES','qr','ML'}
        options.discretization.Sv_uy = Sv_uy;
        options.discretization.Su_vx = Su_vx;
        
%         options.discretization.Diffu_u = Diffu_u;
%         options.discretization.Diffu_v = Diffu_v;       
%         options.discretization.Diffv_u = Diffv_u;
%         options.discretization.Diffv_v = Diffv_v;
     
end


if (order4 == 0)
    options.discretization.Su_vx_BC_lr = Su_vx_BC_lr;
    options.discretization.Su_vx_BC_lu = Su_vx_BC_lu;
    options.discretization.Sv_uy_BC_lr = Sv_uy_BC_lr;
    options.discretization.Sv_uy_BC_lu = Sv_uy_BC_lu;
end

if (order4 == 1)
    options.discretization.Cux3 = Cux3;
    options.discretization.Cuy3 = Cuy3;
    options.discretization.Cvx3 = Cvx3;
    options.discretization.Cvy3 = Cvy3;
    options.discretization.Su_ux_BC3 = Su_ux_BC3;
    options.discretization.Su_uy_BC3 = Su_uy_BC3;
    options.discretization.Sv_vx_BC3 = Sv_vx_BC3;
    options.discretization.Sv_vy_BC3 = Sv_vy_BC3;
    options.discretization.Diffux_div = Diffux_div;
    options.discretization.Diffuy_div = Diffuy_div;
    options.discretization.Diffvx_div = Diffvx_div;
    options.discretization.Diffvy_div = Diffvy_div;
end


%% additional for implicit time stepping diffusion
if (options.time.method==2 && strcmp(visc,'laminar'))
    fprintf(options.output.fcw,'implicit time-stepping for diffusion\n');
    theta = options.time.theta;
    dt    = options.time.dt;
    Omu_inv = options.grid.Omu_inv;
    Omv_inv = options.grid.Omv_inv;
    % implicit time-stepping for diffusion
    % solving (I-dt*Diffu)*uh* = ...
    
    Diffu_impl = speye(Nu,Nu) - theta*dt*spdiags(Omu_inv,0,Nu,Nu)*Diffu;
    Diffv_impl = speye(Nv,Nv) - theta*dt*spdiags(Omv_inv,0,Nv,Nv)*Diffv;
    
    if (options.solversettings.poisson_diffusion == 1)
        % LU decomposition
        fprintf(options.output.fcw,'LU decomposition of diffusion matrices...\n');
        tic;
        [L_diffu, U_diffu] = lu(Diffu_impl);
        [L_diffv, U_diffv] = lu(Diffv_impl);
        toc;
        options.discretization.L_diffu = L_diffu;
        options.discretization.U_diffu = U_diffu;
        options.discretization.L_diffv = L_diffv;
        options.discretization.U_diffv = U_diffv;
    elseif (options.solversettings.poisson_diffusion ==3)
        if (exist(['cg.' mexext],'file')==3)
            [L_diffu, d_diffu] = spdiags(Diffu_impl);
            ndia_diffu  = (length(d_diffu)+1)/2; % number of diagonals needed in cg
            dia_diffu   = d_diffu(ndia_diffu:end);
            L_diffu     = L_diffu(:,ndia_diffu:-1:1);
            
            [L_diffv, d_diffv] = spdiags(Diffv_impl);
            ndia_diffv  = (length(d_diffv)+1)/2; % number of diagonals needed in cg
            dia_diffv   = d_diffv(ndia_diffv:end);
            L_diffv     = L_diffv(:,ndia_diffv:-1:1);
            
        else
            fprintf(options.output.fcw,'No correct CG mex file available, switching to Matlab implementation\n');
            %             poisson = 4;
        end
    end
    
    % CG
    %     [Diffu_diag, d_diffu]   = spdiags(speye(Nu,Nu)-theta*dt*spdiags(Omu_inv,0,Nu,Nu)*Diffu);
    %     nd_diffu                = (length(d_diffu)+1)/2;
    %     dia_diffu               = d_diffu(nd_diffu:end);
    %     B_diffu                 = zeros(Nux_in*Nuy_in,nd_diffu);
    %     B_diffu(:,1:nd_diffu)   = Diffu_diag(:,nd_diffu:-1:1);
    %
    %     [Diffv_diag, d_diffv]   = spdiags(speye(Nv,Nv)-theta*dt*spdiags(Omv_inv,0,Nv,Nv)*Diffv);
    %     nd_diffv                = (length(d_diffv)+1)/2;
    %     dia_diffv               = d_diffv(nd_diffv:end);
    %     B_diffv                 = zeros(Nvx_in*Nvy_in,nd_diffv);
    %     B_diffv(:,1:nd_diffv)   = Diffv_diag(:,nd_diffv:-1:1);
end