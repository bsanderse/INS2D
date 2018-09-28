if (strcmp(visc,'laminar'))
    
    d2u     = Diffu*uh + yDiffu;
    d2v     = Diffv*vh + yDiffv;
    
else
    
%     Diffu_11   = Dux*( nu_ux * 2 * Su_ux) + Duy*( nu_uy * Su_uy);
%     Diffu_12   = Duy*( nu_uy * Sv_uy);
%     Diffv_21   = Dvx*( nu_vx * Su_vx);
%     Diffv_22   = Dvx*( nu_vx * Sv_vx) + Dvy*( nu_vy * 2 * Sv_vy);

    strain_tensor;
    
    C_S  = 0.17;
    filter_length = deltax; % =sqrt(FV size) for uniform grids
    
    % Smagorinsky
%     nu_t = (C_S^2)*(filter_length^2)*sqrt(4*q);
    
    % q-r
    C_d  = deltax^2/8;
    nu_t = C_d * 0.5 * sqrt(4*q) * (1 - alfa / C_d)^2;
    
    nu_t = reshape(nu_t,Npx,Npy);
    
    % ux and vy positions are correct
    % add extra points for periodic BC
    nu_t_ux               = zeros(Nux_in+1,Nuy_in);
    nu_t_ux(2:Nux_in+1,:) = nu_t;
    nu_t_ux(1,:)          = nu_t_ux(Nux_in,:);
    
    nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
    nu_t_vy(:,2:Nvy_in+1) = nu_t;
    nu_t_vy(:,1)          = nu_t_vy(:,Nvy_in);
    
    % construct a larger matrix which can be used for interpolation
    nu_t_ghost            = zeros(Npx+2,Npy+2);
    nu_t_ghost(2:Npx+1,2:Npy+1) = nu_t;
    nu_t_ghost(2:Npx+1,1)   = nu_t(:,end);
    nu_t_ghost(2:Npx+1,end) = nu_t(:,1);
    nu_t_ghost(1,2:Npy+1)   = nu_t(end,:);
    nu_t_ghost(end,2:Npy+1) = nu_t(1,:);
    nu_t_ghost(1,1)         = nu_t(end,end);
    nu_t_ghost(end,end)     = nu_t(1,1);
    nu_t_ghost(1,end)       = nu_t(end,1);
    nu_t_ghost(end,1)       = nu_t(1,end);
    
    xp_ext                  = [xp(1)-hx(1);xp;xp(end)+hx(end)];
    yp_ext                  = [yp(1)-hy(1);yp;yp(end)+hy(end)];
    nu_t_uy                 = interp2(yp_ext',xp_ext,nu_t_ghost,y',xin);
    nu_t_vx                 = interp2(yp_ext',xp_ext,nu_t_ghost,yin',x);

    
    d2u = Dux*((2*(nu + nu_t_ux(:))).*S11(:)) + Duy*((2*(nu + nu_t_uy(:))).*S12(:));
    d2v = Dvx*((2*(nu + nu_t_vx(:))).*S21(:)) + Dvy*((2*(nu + nu_t_vy(:))).*S22(:));

end

