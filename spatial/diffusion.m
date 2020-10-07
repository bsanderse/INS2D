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
        
        
        % compute turbulent viscosity based on S_abs
        switch visc
            case 'LES' % Smagorinsky
                
                C_S  = options.visc.Cs;
                filter_length = deltax; % =sqrt(FV size) for uniform grids
                
                nu_t = (C_S^2)*(filter_length^2)*S_abs;
                
            case 'qr'  % q-r
                % q-r
                C_d  = deltax^2/8;
                nu_t = C_d * 0.5 * S_abs * (1 - alfa / C_d)^2;
                
            case 'ML' % mixing-length
                lm   = options.visc.lm; % mixing length
                nu_t = (lm^2)*S_abs;
                
            otherwise
                error('wrong value for visc parameter');
                
        end
        
        % we now have the turbulent viscosity at all pressure points

        % to compute the diffusion, we need nu_t at ux, uy, vx and vy
        % locations
        % this means we have to reverse the process of strain_tensor.m: go
        % from pressure points back to the ux, uy, vx, vy locations
        [nu_t_ux,nu_t_uy,nu_t_vx,nu_t_vy] = interpolate_nu(nu_t,options);
        
        % now the total diffusive terms (laminar + turbulent) is as follows
        % note that the factor 2 is because
        % tau = 2*(nu+nu_t)*S(u), with S(u) = 0.5*(grad u + (grad u)^T)
        Dux = options.discretization.Dux;
        Duy = options.discretization.Duy;
        Dvx = options.discretization.Dvx;
        Dvy = options.discretization.Dvy;
        
        nu  = 1/options.fluid.Re; % molecular viscosity
        
        d2u = Dux*(2*(nu + nu_t_ux).*S11(:)) + Duy*(2*(nu + nu_t_uy).*S12(:));
        d2v = Dvx*(2*(nu + nu_t_vx).*S21(:)) + Dvy*(2*(nu + nu_t_vy).*S22(:));
        
        
        if (getJacobian == 1)
            % freeze nu_t, i.e. we skip the derivative of nu_t wrt V in 
            % the Jacobian
            Su_ux = options.discretization.Su_ux;
            Su_uy = options.discretization.Su_uy;
            Su_vx = options.discretization.Su_vx;
            Sv_vx = options.discretization.Sv_vx;
            Sv_vy = options.discretization.Sv_vy;
            Sv_uy = options.discretization.Sv_uy;
            
            % if (getJacobian == 0)
            %     Jacu = spalloc(Nu,Nu+Nv,0);
            %     Jacv = spalloc(Nv,Nu+Nv,0);    
            % elseif (getJacobian == 1)
            %     Jacu = [Su_ux + (1/2)*Su_uy        (1/2)*Sv_uy];
            %     Jacv = [(1/2)*Su_vx          (1/2)*Sv_vx + Sv_vy];
            % end
            
            N1 = options.grid.N1;
            N2 = options.grid.N2;
            N3 = options.grid.N3;
            N4 = options.grid.N4;
            Jacu    = [Dux*(2*spdiags(nu + nu_t_ux,0,N1,N1))*Su_ux + ...
                          Duy*(2*spdiags(nu + nu_t_uy,0,N2,N2))*(1/2)*Su_uy ... 
                       Duy*(2*spdiags(nu + nu_t_uy,0,N2,N2))*(1/2)*Sv_uy];
            Jacv    = [Dvx*(2*spdiags(nu + nu_t_vx,0,N3,N3))*(1/2)*Su_vx ...
                       Dvx*(2*spdiags(nu + nu_t_vx,0,N3,N3))*(1/2)*Sv_vx + ...
                          Dvy*(2*spdiags(nu + nu_t_vy,0,N4,N4))*Sv_vy];
        end
                
    case 'keps' % (k-e)
        
        error('k-e implementation in diffusion.m not finished');
        
    otherwise
        
        error('wrong specification of viscosity model');
        
end


end

