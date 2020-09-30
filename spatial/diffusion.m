function [d2u, d2v, Jacu, Jacv] = diffusion(V,t,options,getJacobian)
% evaluate diffusive terms and optionally Jacobian

visc = options.case.visc;


Nu = options.grid.Nu;
Nv = options.grid.Nv;

indu = options.grid.indu;
indv = options.grid.indv;
uh   = V(indu);
vh   = V(indv);



Jacu = spalloc(Nu,Nu+Nv,0);
Jacv = spalloc(Nv,Nu+Nv,0);

switch visc
    
    case 'laminar'
        
        Diffu  = options.discretization.Diffu;
        Diffv  = options.discretization.Diffv;
        yDiffu = options.discretization.yDiffu;
        yDiffv = options.discretization.yDiffv;
        
        d2u     = Diffu*uh + yDiffu;
        d2v     = Diffv*vh + yDiffv;
        
        if (getJacobian == 1)
            Jacu    = [Diffu spalloc(Nu,Nv,0)];
            Jacv    = [spalloc(Nv,Nu,0) Diffv];
        end
        
    case {'qr','LES','ML'}
        
        % get components of strain tensor and its magnitude;
        % the magnitude S_abs is evaluated at pressure points
        [S11,S12,S21,S22,S_abs] = strain_tensor(V,t,options,getJacobian);
        
        Dux = options.discretization.Dux;
        Duy = options.discretization.Duy;
        Dvx = options.discretization.Dvx;
        Dvy = options.discretization.Dvy;
        
        %
        switch visc
            case 'LES' % Smagorinsky
                
                C_S  = 0.17;
                filter_length = deltax; % =sqrt(FV size) for uniform grids
                
                nu_t = (C_S^2)*(filter_length^2)*S_abs;
                
            case 'qr'  % q-r
                % q-r
                C_d  = deltax^2/8;
                nu_t = C_d * 0.5 * S_abs * (1 - alfa / C_d)^2;
                
            case 'ML' % mixing-length
                lm   = 10; % mixing length
                nu_t = (lm^2)*S_abs;
        end
        
        % we now have the turbulent viscosity at all pressure points
        Npx = options.grid.Npx;
        Npy = options.grid.Npy;
        Nux_in = options.grid.Nux_in;
        Nuy_in = options.grid.Nuy_in;
        Nvx_in = options.grid.Nvx_in;
        Nvy_in = options.grid.Nvy_in;        
        
        nu_t = reshape(nu_t,Npx,Npy);
        
        % to compute the diffusion, we need nu_t at ux, uy, vx and vy
        % locations
        
        % this means we have to reverse the process of strain_tensor.m: go
        % from pressure points back to the ux, uy, vx, vy locations
        
        BC = options.BC;
        if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
            % ux positions are correct
            % add extra points for periodic BC
            nu_t_ux               = zeros(Nux_in+1,Nuy_in);
            nu_t_ux(2:Nux_in+1,:) = nu_t;
            nu_t_ux(1,:)          = nu_t_ux(Nux_in,:); % bit unclear why 1 and not Nux_in+1, see also strain_tensor.m
            
        elseif (strcmp(BC.u.left,'dir') && strcmp(BC.u.right,'pres'))
            % ux positions are correct
            % add extra points for outflow BC
            nu_t_ux               = zeros(Nux_in+1,Nuy_in);
            nu_t_ux(1:Nux_in,:)   = nu_t;
            nu_t_ux(Nux_in+1,:)   = nu_t_ux(Nux_in,:);            
            
        else
            error('BC not implemented in strain_tensor.m');
        end
        
        
        if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
            % vy positions are correct            
            nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
            nu_t_vy(:,2:Nvy_in+1) = nu_t;
            nu_t_vy(:,1)          = nu_t_vy(:,Nvy_in); % bit unclear why 1 and not Nvy_in+1, see also strain_tensor.m
            
        elseif (strcmp(BC.v.low,'pres') && strcmp(BC.v.up,'pres'))
            % vy positions are correct    
            % add extra points for outflow BC
            
            nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
            nu_t_vy(:,2:Nvy_in)   = nu_t;
            nu_t_vy(:,1)          = nu_t_vy(:,2); 
            nu_t_vy(:,Nvy_in+1)   = nu_t_vy(:,Nvy_in);

            

            
        else
            error('BC not implemented in strain_tensor.m');
        end
        
        
        % for uy and vx, the locations are simply given by (xin,y) and
        % (x,yin)
        % we therefore create nu_t_extended:
%         nu_t_uy = interp2(yp',xp,nu_t,y',xin);
%         nu_t_vx = interp2(yp',xp,nu_t,yin',x);
%         
%         
        
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
        %
        
        
        %     nu_t = reshape(nu_t,Npx,Npy);
        %         S11 = (1/2)* 2*(Su_ux*uh + ySu_ux);
        %         S12 = (1/2)* (Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy);
        %         S21 = (1/2)* (Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx);
        %         S22 = (1/2)* 2*(Sv_vy*vh + ySv_vy);
        d2u = Dux*((2*(nu + nu_t_ux(:))).*S11(:)) + Duy*((2*(nu + nu_t_uy(:))).*S12(:));
        d2v = Dvx*((2*(nu + nu_t_vx(:))).*S21(:)) + Dvy*((2*(nu + nu_t_vy(:))).*S22(:));
        
        
        % code snippet from k-eps equations:
        %       nu_ux      = spdiags( nu + Cmu* kAux.^2 ./ eAux,0,N1,N1);
        %       nu_uy      = spdiags( nu + Cmu* kAuy.^2 ./ eAuy,0,N2,N2);
        %       nu_vx      = spdiags( nu + Cmu* kAvx.^2 ./ eAvx,0,N3,N3);
        %       nu_vy      = spdiags( nu + Cmu* kAvy.^2 ./ eAvy,0,N4,N4);
        %
        %       % difference u and v to ux, uy, vx, vy positions
        %       uSux       = Su_ux*uh + ySu_ux;
        %       uSuy       = Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy;
        %       vSvx       = Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx;
        %       vSvy       = Sv_vy*vh + ySv_vy;
        %
        %       % residual
        %       yDiffu     = Dux*( nu_ux* 2 * uSux ) + Duy*( nu_uy* uSuy );
        %       yDiffv     = Dvx*( nu_vx* vSvx ) + Dvy*( nu_vy* 2 * vSvy);
        
        error('implementation in diffusion.m not finished');

        
    case 'turbulent' % (k-e)
        
        error('k-e implementation in diffusion.m not finished');
        
    otherwise
        
        error('wrong specification of viscosity model');
        
end

% else % 'turbulent' (k-e) or 'LES'
%
%     %     Diffu_11   = Dux*( nu_ux * 2 * Su_ux) + Duy*( nu_uy * Su_uy);
%     %     Diffu_12   = Duy*( nu_uy * Sv_uy);
%     %     Diffv_21   = Dvx*( nu_vx * Su_vx);
%     %     Diffv_22   = Dvx*( nu_vx * Sv_vx) + Dvy*( nu_vy * 2 * Sv_vy);
%
%     strain_tensor;
%
%     C_S  = 0.17;
%     filter_length = deltax; % =sqrt(FV size) for uniform grids
%
%     % Smagorinsky
%     %     nu_t = (C_S^2)*(filter_length^2)*sqrt(4*q);
%
%     % q-r
%     C_d  = deltax^2/8;
%     nu_t = C_d * 0.5 * sqrt(4*q) * (1 - alfa / C_d)^2;
%
%     nu_t = reshape(nu_t,Npx,Npy);
%
%     % ux and vy positions are correct
%     % add extra points for periodic BC
%     nu_t_ux               = zeros(Nux_in+1,Nuy_in);
%     nu_t_ux(2:Nux_in+1,:) = nu_t;
%     nu_t_ux(1,:)          = nu_t_ux(Nux_in,:);
%
%     nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
%     nu_t_vy(:,2:Nvy_in+1) = nu_t;
%     nu_t_vy(:,1)          = nu_t_vy(:,Nvy_in);
%
%     % construct a larger matrix which can be used for interpolation
%     nu_t_ghost            = zeros(Npx+2,Npy+2);
%     nu_t_ghost(2:Npx+1,2:Npy+1) = nu_t;
%     nu_t_ghost(2:Npx+1,1)   = nu_t(:,end);
%     nu_t_ghost(2:Npx+1,end) = nu_t(:,1);
%     nu_t_ghost(1,2:Npy+1)   = nu_t(end,:);
%     nu_t_ghost(end,2:Npy+1) = nu_t(1,:);
%     nu_t_ghost(1,1)         = nu_t(end,end);
%     nu_t_ghost(end,end)     = nu_t(1,1);
%     nu_t_ghost(1,end)       = nu_t(end,1);
%     nu_t_ghost(end,1)       = nu_t(1,end);
%
%     xp_ext                  = [xp(1)-hx(1);xp;xp(end)+hx(end)];
%     yp_ext                  = [yp(1)-hy(1);yp;yp(end)+hy(end)];
%     nu_t_uy                 = interp2(yp_ext',xp_ext,nu_t_ghost,y',xin);
%     nu_t_vx                 = interp2(yp_ext',xp_ext,nu_t_ghost,yin',x);
%
%
%     d2u = Dux*((2*(nu + nu_t_ux(:))).*S11(:)) + Duy*((2*(nu + nu_t_uy(:))).*S12(:));
%     d2v = Dvx*((2*(nu + nu_t_vx(:))).*S21(:)) + Dvy*((2*(nu + nu_t_vy(:))).*S22(:));
%
% end

end

