function options = set_bc_vectors(t,options)
% construct boundary conditions


%% get settings from options structure

% steady
steady = options.case.steady;
% 4th order
order4 = options.discretization.order4;
% Reynolds number
Re = options.fluid.Re;
% boundary conditions
BC = options.BC;

global uBC vBC dudtBC dvdtBC;

% type of stress tensor
visc = options.case.visc;


if (order4 == 1)
    alfa   = options.discretization.alfa;
end

%% grid settings

% number of interior points and boundary points
Nux_in = options.grid.Nux_in;
Nvy_in = options.grid.Nvy_in;
Np     = options.grid.Np;

xin = options.grid.xin;
yin = options.grid.yin;
x   = options.grid.x;
y   = options.grid.y;
hx  = options.grid.hx;
hy  = options.grid.hy;
xp  = options.grid.xp;
yp  = options.grid.yp;


%% get BC values
uLo      = uBC(x,y(1),t,options);
uUp      = uBC(x,y(end),t,options);
% uLe      = uBC(x(1),y,t,options);
% uRi      = uBC(x(end),y,t,options);

uLo_i    = uBC(xin,y(1),t,options);
uUp_i    = uBC(xin,y(end),t,options);
uLe_i    = uBC(x(1),yp,t,options);
uRi_i    = uBC(x(end),yp,t,options);

% vLo      = vBC(x,y(1),t,options);
% vUp      = vBC(x,y(end),t,options);
vLe      = vBC(x(1),y,t,options);
vRi      = vBC(x(end),y,t,options);

vLo_i    = vBC(xp,y(1),t,options);
vUp_i    = vBC(xp,y(end),t,options);
vLe_i    = vBC(x(1),yin,t,options);
vRi_i    = vBC(x(end),yin,t,options);

if (steady == 0 && options.BC.BC_unsteady==1)
    dudtLe_i   = dudtBC(x(1),options.grid.yp,t,options);
    dudtRi_i   = dudtBC(x(end),options.grid.yp,t,options);
    
    dvdtLo_i   = dvdtBC(options.grid.xp,y(1),t,options);
    dvdtUp_i   = dvdtBC(options.grid.xp,y(end),t,options);
end

pLe = options.BC.pLe;
pRi = options.BC.pRi;
pLo = options.BC.pLo;
pUp = options.BC.pUp;


%% boundary conditions for divergence

Mx_BC = options.discretization.Mx_BC;
My_BC = options.discretization.My_BC;


if (order4==1)
    Mx_BC3 = options.discretization.Mx_BC3;
    My_BC3 = options.discretization.My_BC3;
end


% Mx
ybc    = kron(uLe_i,Mx_BC.ybc1) + kron(uRi_i,Mx_BC.ybc2);
yMx    = Mx_BC.Bbc*ybc;
if (order4==1)
    ybc3 = kron(uLe_i,Mx_BC3.ybc1) + kron(uRi_i,Mx_BC3.ybc2);
    yMx3 = Mx_BC3.Bbc*ybc3;
    yMx  = alfa*yMx - yMx3;
end

% My
ybc    = kron(My_BC.ybc1,vLo_i) + kron(My_BC.ybc2,vUp_i);
yMy    = My_BC.Bbc*ybc;
if (order4==1)
    ybc3 = kron(My_BC3.ybc1,vLo_i) + kron(My_BC3.ybc2,vUp_i);
    yMy3 = My_BC3.Bbc*ybc3;
    yMy  = alfa*yMy - yMy3;
end

yM     = yMx + yMy;
options.discretization.yM  = yM;


%% time derivative of divergence
if (steady == 0)
    if (options.BC.BC_unsteady==1)
        ybc    = kron(dudtLe_i,Mx_BC.ybc1) + kron(dudtRi_i,Mx_BC.ybc2);
        ydMx   = Mx_BC.Bbc*ybc;
        if (order4==1)
            ybc3  = kron(dudtLe_i,Mx_BC3.ybc1) + kron(dudtRi_i,Mx_BC3.ybc2);
            ydMx3 = Mx_BC3.Bbc*ybc3;
            ydMx  = alfa*ydMx - ydMx3;
        end
        
        % My
        ybc    = kron(My_BC.ybc1,dvdtLo_i) + kron(My_BC.ybc2,dvdtUp_i);
        ydMy   = My_BC.Bbc*ybc;
        if (order4==1)
            ybc3  = kron(My_BC3.ybc1,dvdtLo_i) + kron(My_BC3.ybc2,dvdtUp_i);
            ydMy3 = My_BC3.Bbc*ybc3;
            ydMy  = alfa*ydMy - ydMy3;
        end
        
        ydM    = ydMx + ydMy;
        
        options.discretization.ydM = ydM;
    else
        options.discretization.ydM = zeros(Np,1);
    end
end

% if (ibm==1)
%     ydM = [ydM; zeros(n_ibm,1)];
%     yM = [yM; zeros(n_ibm,1)];
% end

%% boundary conditions for pressure

% left and right side
y1D_le            = zeros(Nux_in,1);
y1D_ri            = zeros(Nux_in,1);
if (strcmp(BC.u.left,'pres'))
    y1D_le(1)     = -1;
end
if (strcmp(BC.u.right,'pres'))
    y1D_ri(end)   = 1;
end
y_px              = kron(hy.*pLe,y1D_le) + kron(hy.*pRi,y1D_ri);

% lower and upper side
y1D_lo            = zeros(Nvy_in,1);
y1D_up            = zeros(Nvy_in,1);
if (strcmp(BC.v.low,'pres'))
    y1D_lo(1)     = -1;
end
if (strcmp(BC.v.up,'pres'))
    y1D_up(end)   = 1;
end
y_py              = kron(y1D_lo,hx.*pLo) + kron(y1D_up,hx.*pUp);


options.discretization.y_px = y_px;
options.discretization.y_py = y_py;


%% boundary conditions for averaging

Au_ux_BC = options.discretization.Au_ux_BC;
Au_uy_BC = options.discretization.Au_uy_BC;
Av_vx_BC = options.discretization.Av_vx_BC;
Av_vy_BC = options.discretization.Av_vy_BC;
if (order4==1)
    Au_ux_BC3 = options.discretization.Au_ux_BC3;
    Au_uy_BC3 = options.discretization.Au_uy_BC3;
    Av_vx_BC3 = options.discretization.Av_vx_BC3;
    Av_vy_BC3 = options.discretization.Av_vy_BC3;
end

% Au_ux
% uLe_i  = interp1(y,uLe,yp);
% uRi_i  = interp1(y,uRi,yp);
ybc     = kron(uLe_i,Au_ux_BC.ybc1) + kron(uRi_i,Au_ux_BC.ybc2);
yAu_ux  = Au_ux_BC.Bbc*ybc;
if (order4==1)
    ybc3    = kron(uLe_i,Au_ux_BC3.ybc1) + kron(uRi_i,Au_ux_BC3.ybc2);
    yAu_ux3 = Au_ux_BC3.Bbc*ybc3;
end

% Au_uy
ybc     = kron(Au_uy_BC.ybc1,uLo_i) + kron(Au_uy_BC.ybc2,uUp_i);
yAu_uy  = Au_uy_BC.Bbc*ybc;
if (order4==1)
    ybc3    = kron(Au_uy_BC3.ybc1,uLo_i) + kron(Au_uy_BC3.ybc2,uUp_i);
    yAu_uy3 = Au_uy_BC3.Bbc*ybc3;
end

% Av_vx
ybc     = kron(vLe_i,Av_vx_BC.ybc1) + kron(vRi_i,Av_vx_BC.ybc2);
yAv_vx  = Av_vx_BC.Bbc*ybc;
if (order4==1)
    ybc3    = kron(vLe_i,Av_vx_BC3.ybc1) + kron(vRi_i,Av_vx_BC3.ybc2);
    yAv_vx3 = Av_vx_BC3.Bbc*ybc3;
end

% Av_vy
% vLo_i  = interp1(x,vLo,xp);
% vUp_i  = interp1(x,vUp,xp);
ybc    = kron(Av_vy_BC.ybc1,vLo_i) + kron(Av_vy_BC.ybc2,vUp_i);
yAv_vy = Av_vy_BC.Bbc*ybc;
if (order4==1)
    ybc3    = kron(Av_vy_BC3.ybc1,vLo_i) + kron(Av_vy_BC3.ybc2,vUp_i);
    yAv_vy3 = Av_vy_BC3.Bbc*ybc3;
end

options.discretization.yAu_ux = yAu_ux;
options.discretization.yAu_uy = yAu_uy;
options.discretization.yAv_vx = yAv_vx;
options.discretization.yAv_vy = yAv_vy;

if (order4==1)
    options.discretization.yAu_ux3 = yAu_ux3;
    options.discretization.yAu_uy3 = yAu_uy3;
    options.discretization.yAv_vx3 = yAv_vx3;
    options.discretization.yAv_vy3 = yAv_vy3;
end



%% boundary conditions for diffusion

Su_ux_BC = options.discretization.Su_ux_BC;
Su_uy_BC = options.discretization.Su_uy_BC;
Sv_vx_BC = options.discretization.Sv_vx_BC;
Sv_vy_BC = options.discretization.Sv_vy_BC;
Dux = options.discretization.Dux;
Duy = options.discretization.Duy;
Dvx = options.discretization.Dvx;
Dvy = options.discretization.Dvy;


if (order4==0)
    
    Su_vx_BC_lr = options.discretization.Su_vx_BC_lr;
    Su_vx_BC_lu = options.discretization.Su_vx_BC_lu;
    Sv_uy_BC_lr = options.discretization.Sv_uy_BC_lr;
    Sv_uy_BC_lu = options.discretization.Sv_uy_BC_lu;
    
    % Su_ux
    % uLe_i  = interp1(y,uLe,yp);
    % uRi_i  = interp1(y,uRi,yp);
    ybc    = kron(uLe_i,Su_ux_BC.ybc1) + kron(uRi_i,Su_ux_BC.ybc2);
    ySu_ux = Su_ux_BC.Bbc*ybc;
    
    % Su_uy
    % uLo_i  = interp1(x,uLo,xin);
    % uUp_i  = interp1(x,uUp,xin);
    ybc    = kron(Su_uy_BC.ybc1,uLo_i) + kron(Su_uy_BC.ybc2,uUp_i);
    ySu_uy = Su_uy_BC.Bbc*ybc;
    
    % Sv_uy
    % left/right
    ybc       = kron(vLe,Sv_uy_BC_lr.ybc1) + kron(vRi,Sv_uy_BC_lr.ybc2);
    ySv_uy_lr = Sv_uy_BC_lr.Bbc*ybc;
    % low/up
    % vLo_i     = interp1(x,vLo,xp);
    % vUp_i     = interp1(x,vUp,xp);
    ybc       = kron(Sv_uy_BC_lu.ybc1,vLo_i) + kron(Sv_uy_BC_lu.ybc2,vUp_i);
    ySv_uy_lu = Sv_uy_BC_lr.B2D*Sv_uy_BC_lu.Bbc*ybc;
    
    ySv_uy    = ySv_uy_lr + ySv_uy_lu;
    
    % Su_vx
    % low/up
    ybc       = kron(Su_vx_BC_lu.ybc1,uLo) + kron(Su_vx_BC_lu.ybc2,uUp);
    ySu_vx_lu = Su_vx_BC_lu.Bbc*ybc;
    % left/right
    % uLe_i     = interp1(y,uLe,yp);
    % uRi_i     = interp1(y,uRi,yp);
    ybc       = kron(uLe_i,Su_vx_BC_lr.ybc1) + kron(uRi_i,Su_vx_BC_lr.ybc2);
    ySu_vx_lr = Su_vx_BC_lu.B2D*Su_vx_BC_lr.Bbc*ybc;
    ySu_vx    = ySu_vx_lr + ySu_vx_lu;
    
    % Sv_vx
    % vLe_i  = interp1(y,vLe,yin);
    % vRi_i  = interp1(y,vRi,yin);
    ybc    = kron(vLe_i,Sv_vx_BC.ybc1) + kron(vRi_i,Sv_vx_BC.ybc2);
    ySv_vx = Sv_vx_BC.Bbc*ybc;
    
    % Sv_vy
    % vLo_i  = interp1(x,vLo,xp);
    % vUp_i  = interp1(x,vUp,xp);
    ybc    = kron(Sv_vy_BC.ybc1,vLo_i) + kron(Sv_vy_BC.ybc2,vUp_i);
    ySv_vy = Sv_vy_BC.Bbc*ybc;
    
    
    if ( strcmp(visc,'laminar') )
        
        yDiffu = Dux*( (1/Re)* ySu_ux) + Duy*( (1/Re)* ySu_uy);
        yDiffv = Dvx*( (1/Re)* ySv_vx) + Dvy*( (1/Re)* ySv_vy);
        
    elseif ( strcmp(visc,'turbulent') )
        
        yDiffu  = Dux*( (1/Re) * 2*ySu_ux) + Duy*( (1/Re) * ySu_uy + (1/Re) * ySv_uy);
        yDiffv  = Dvx*( (1/Re) * ySv_vx + (1/Re) * ySu_vx) + Dvy*( (1/Re) * 2*ySv_vy);
        
    end
    
end

if (order4==1)
    
    Su_ux_BC3 = options.discretization.Su_ux_BC3;
    Su_uy_BC3 = options.discretization.Su_uy_BC3;
    Sv_vx_BC3 = options.discretization.Sv_vx_BC3;
    Sv_vy_BC3 = options.discretization.Sv_vy_BC3;
    Diffux_div = options.discretization.Diffux_div;
    Diffuy_div = options.discretization.Diffuy_div;
    Diffvx_div = options.discretization.Diffvx_div;
    Diffvy_div = options.discretization.Diffvy_div;
    
    ybc1    = kron(uLe_i,Su_ux_BC.ybc1) + kron(uRi_i,Su_ux_BC.ybc2);
    ybc3    = kron(uLe_i,Su_ux_BC3.ybc1) + kron(uRi_i,Su_ux_BC3.ybc2);
    ySu_ux  = alfa*Su_ux_BC.Bbc*ybc1 - Su_ux_BC3.Bbc*ybc3;
    
    ybc1    = kron(Su_uy_BC.ybc1,uLo_i) + kron(Su_uy_BC.ybc2,uUp_i);
    ybc3    = kron(Su_uy_BC3.ybc1,uLo_i) + kron(Su_uy_BC3.ybc2,uUp_i);
    ySu_uy  = alfa*Su_uy_BC.Bbc*ybc1 - Su_uy_BC3.Bbc*ybc3;
    
    ybc1    = kron(vLe_i,Sv_vx_BC.ybc1) + kron(vRi_i,Sv_vx_BC.ybc2);
    ybc3    = kron(vLe_i,Sv_vx_BC3.ybc1) + kron(vRi_i,Sv_vx_BC3.ybc2);
    ySv_vx  = alfa*Sv_vx_BC.Bbc*ybc1 - Sv_vx_BC3.Bbc*ybc3;
    
    ybc1    = kron(Sv_vy_BC.ybc1,vLo_i) + kron(Sv_vy_BC.ybc2,vUp_i);
    ybc3    = kron(Sv_vy_BC3.ybc1,vLo_i) + kron(Sv_vy_BC3.ybc2,vUp_i);
    ySv_vy  = alfa*Sv_vy_BC.Bbc*ybc1 - Sv_vy_BC3.Bbc*ybc3;
    
    yDiffu  = (1/Re)*(Diffux_div*ySu_ux + Diffuy_div*ySu_uy);
    yDiffv  = (1/Re)*(Diffvx_div*ySv_vx + Diffvy_div*ySv_vy);
end

options.discretization.yDiffu = yDiffu;
options.discretization.yDiffv = yDiffv;


%% boundary conditions for interpolation

Iu_ux_BC    = options.discretization.Iu_ux_BC;
Iv_uy_BC_lr = options.discretization.Iv_uy_BC_lr;
Iv_uy_BC_lu = options.discretization.Iv_uy_BC_lu;
Iu_vx_BC_lr = options.discretization.Iu_vx_BC_lr;
Iu_vx_BC_lu = options.discretization.Iu_vx_BC_lu;
Iv_vy_BC    = options.discretization.Iv_vy_BC;

if (order4 == 1)
    Iu_ux_BC3    = options.discretization.Iu_ux_BC3;
    Iv_uy_BC_lu3 = options.discretization.Iv_uy_BC_lu3;
    Iv_uy_BC_lr3 = options.discretization.Iv_uy_BC_lr3;
    Iu_vx_BC_lu3 = options.discretization.Iu_vx_BC_lu3;
    Iu_vx_BC_lr3 = options.discretization.Iu_vx_BC_lr3;
    Iv_vy_BC3    = options.discretization.Iv_vy_BC3;
end

% Iu_ux
% uLe_i           = interp1(y,uLe,yp);
% uRi_i           = interp1(y,uRi,yp);
ybc             = kron(uLe_i,Iu_ux_BC.ybc1) + kron(uRi_i,Iu_ux_BC.ybc2);
yIu_ux          = Iu_ux_BC.Bbc*ybc;
if (order4==1)
    ybc3         = kron(uLe_i,Iu_ux_BC3.ybc1) + kron(uRi_i,Iu_ux_BC3.ybc2);
    yIu_ux3      = Iu_ux_BC3.Bbc*ybc3;
end

% Iv_uy
% left/right
ybc             = kron(vLe,Iv_uy_BC_lr.ybc1) + kron(vRi,Iv_uy_BC_lr.ybc2);
yIv_uy_lr       = Iv_uy_BC_lr.Bbc*ybc;
% low/up
% vLo_i           = interp1(x,vLo,xp);
% vUp_i           = interp1(x,vUp,xp);
ybc             = kron(Iv_uy_BC_lu.ybc1,vLo_i) + kron(Iv_uy_BC_lu.ybc2,vUp_i);
yIv_uy_lu       = Iv_uy_BC_lr.B2D*Iv_uy_BC_lu.Bbc*ybc;
yIv_uy          = yIv_uy_lr + yIv_uy_lu;

if (order4==1)
    if (strcmp(BC.v.low,'dir'))
        vLe_ext     = [2*vLe(1)-vLe(2);vLe];
        vRi_ext     = [2*vRi(1)-vRi(2);vRi];
    elseif (strcmp(BC.v.low,'per'))
        vLe_ext     = [0;vLe];
        vRi_ext     = [0;vRi];
    elseif (strcmp(BC.v.low,'pres'))
        vLe_ext     = [vLe(2);vLe]; % zero gradient
        vRi_ext     = [vRi(2);vRi]; % zero gradient
    end
    if (strcmp(BC.v.up,'dir'))
        vLe_ext     = [vLe_ext;2*vLe(end)-vLe(end-1)];
        vRi_ext     = [vRi_ext;2*vRi(1)-vRi(2)];
    elseif (strcmp(BC.v.up,'per'))
        vLe_ext     = [vLe_ext;0];
        vRi_ext     = [vRi_ext;0];
    elseif (strcmp(BC.v.up,'pres'))
        vLe_ext     = [vLe_ext;vLe(end-1)]; % zero gradient
        vRi_ext     = [vRi_ext;vRi(end-1)]; % zero gradient
    end
    ybc3         = kron(vLe_ext,Iv_uy_BC_lr3.ybc1) + kron(vRi_ext,Iv_uy_BC_lr3.ybc2);
    yIv_uy_lr3   = Iv_uy_BC_lr3.Bbc*ybc3;
    
    ybc3         = kron(Iv_uy_BC_lu3.ybc1,vLo_i) + kron(Iv_uy_BC_lu3.ybc2,vUp_i);
    yIv_uy_lu3   = Iv_uy_BC_lr3.B2D*Iv_uy_BC_lu3.Bbc*ybc3;
    yIv_uy3      = yIv_uy_lr3 + yIv_uy_lu3;
    
end


% Iu_vx
% low/up
ybc             = kron(Iu_vx_BC_lu.ybc1,uLo) + kron(Iu_vx_BC_lu.ybc2,uUp);
yIu_vx_lu       = Iu_vx_BC_lu.Bbc*ybc;
% left/right
% uLe_i           = interp1(y,uLe,yp);
% uRi_i           = interp1(y,uRi,yp);
ybc             = kron(uLe_i,Iu_vx_BC_lr.ybc1) + kron(uRi_i,Iu_vx_BC_lr.ybc2);
yIu_vx_lr       = Iu_vx_BC_lu.B2D*Iu_vx_BC_lr.Bbc*ybc;
yIu_vx          = yIu_vx_lr + yIu_vx_lu;

if (order4==1)
    if (strcmp(BC.u.left,'dir'))
        uLo_ext     = [2*uLo(1)-uLo(2);uLo];
        uUp_ext     = [2*uUp(1)-uUp(2);uUp];
    elseif (strcmp(BC.u.left,'per'))
        uLo_ext     = [0;uLo];
        uUp_ext     = [0;uUp];
    elseif (strcmp(BC.u.left,'pres'))
        uLo_ext     = [uLo(2);uLo]; % zero gradient
        uUp_ext     = [uUp(2);uUp]; % zero gradient
    end
    if (strcmp(BC.u.right,'dir'))
        uLo_ext     = [uLo_ext;2*uLo(end)-uLo(end-1)];
        uUp_ext     = [uUp_ext;2*uUp(1)-uUp(2)];
    elseif (strcmp(BC.u.right,'per'))
        uLo_ext     = [uLo_ext;0];
        uUp_ext     = [uUp_ext;0];
    elseif (strcmp(BC.u.right,'pres'))
        uLo_ext     = [uLo_ext;uLo(end-1)]; % zero gradient
        uUp_ext     = [uUp_ext;uUp(end-1)]; % zero gradient
    end
    ybc3         = kron(Iu_vx_BC_lu3.ybc1,uLo_ext) + kron(Iu_vx_BC_lu3.ybc2,uUp_ext);
    yIu_vx_lu3   = Iu_vx_BC_lu3.Bbc*ybc3;
    
    ybc3         = kron(uLe_i,Iu_vx_BC_lr3.ybc1) + kron(uRi_i,Iu_vx_BC_lr3.ybc2);
    yIu_vx_lr3   = Iu_vx_BC_lu3.B2D*Iu_vx_BC_lr3.Bbc*ybc3;
    yIu_vx3      = yIu_vx_lr3 + yIu_vx_lu3;
    
end

% Iv_vy
% vLo_i           = interp1(x,vLo,xp);
% vUp_i           = interp1(x,vUp,xp);
ybc             = kron(Iv_vy_BC.ybc1,vLo_i) + kron(Iv_vy_BC.ybc2,vUp_i);
yIv_vy          = Iv_vy_BC.Bbc*ybc;
if (order4==1)
    ybc3         = kron(Iv_vy_BC3.ybc1,vLo_i) + kron(Iv_vy_BC3.ybc2,vUp_i);
    yIv_vy3      = Iv_vy_BC3.Bbc*ybc3;
end

options.discretization.yIu_ux = yIu_ux;
options.discretization.yIv_uy = yIv_uy;
options.discretization.yIu_vx = yIu_vx;
options.discretization.yIv_vy = yIv_vy;

if (order4==1)
    options.discretization.yIu_ux3 = yIu_ux3;
    options.discretization.yIv_uy3 = yIv_uy3;
    options.discretization.yIu_vx3 = yIu_vx3;
    options.discretization.yIv_vy3 = yIv_vy3;
end


end