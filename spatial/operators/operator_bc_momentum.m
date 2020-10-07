%% THIS FILE IS OBSOLETE NOW THAT WE USE set_bc_vectors


%% operator for boundary conditions


%% pressure

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


%% averaging
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


%% diffusion

if (order4==0)
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
    
    switch visc
        
        case 'laminar'
        
            yDiffu = Dux*( (1/Re)* ySu_ux) + Duy*( (1/Re)* ySu_uy);
            yDiffv = Dvx*( (1/Re)* ySv_vx) + Dvy*( (1/Re)* ySv_vy);
        
        case 'keps'
        
            yDiffu  = Dux*( (1/Re) * 2*ySu_ux) + Duy*( (1/Re) * ySu_uy + (1/Re) * ySv_uy);
            yDiffv  = Dvx*( (1/Re) * ySv_vx + (1/Re) * ySu_vx) + Dvy*( (1/Re) * 2*ySv_vy);
        
    end
    
end

if (order4==1)
    
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

%% interpolation
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
