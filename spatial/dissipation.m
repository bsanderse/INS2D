function [Phi, Jac] = dissipation(V,t,options,getJacobian)
% evaluate dissipation terms and optionally Jacobian
% the dissipation terms are defined on the pressure grid (where the kinetic
% energy is defined) and are a function of the velocity field V

% note the similarity with strain_tensor.m

visc = options.case.visc;

Nx = options.grid.Nx;
Ny = options.grid.Ny;

Np = options.grid.Np;
Npx = options.grid.Npx;
Npy = options.grid.Npy;


Nu = options.grid.Nu;
Nux_in = options.grid.Nux_in;
Nux_b  = options.grid.Nux_b;
Nux_t  = options.grid.Nux_t;
Nuy_in = options.grid.Nuy_in;
Nuy_b  = options.grid.Nuy_b;
Nuy_t  = options.grid.Nuy_t;

Nv = options.grid.Nv;
Nvx_in = options.grid.Nvx_in;
Nvx_b  = options.grid.Nvx_b;
Nvx_t  = options.grid.Nvx_t;
Nvy_in = options.grid.Nvy_in;
Nvy_b  = options.grid.Nvy_b;
Nvy_t  = options.grid.Nvy_t;

Omu = options.grid.Omu;
Omv = options.grid.Omv;

BC = options.BC;

hx = options.grid.hx;
hy = options.grid.hy;
gx = options.grid.gx;
gy = options.grid.gy;

indu = options.grid.indu;
indv = options.grid.indv;
uh   = V(indu);
vh   = V(indv);

Jac = spalloc(Np,Nu+Nv,0);

switch visc
    
    case 'laminar'
        
        %         Diffu  = options.discretization.Diffu;
        %         Diffv  = options.discretization.Diffv;
        %         yDiffu = options.discretization.yDiffu;
        %         yDiffv = options.discretization.yDiffv;
        %
        %         d2u     = Diffu*uh + yDiffu;
        %         d2v     = Diffv*vh + yDiffv;
        
        % get du/dx at the ux locations
        Su_ux = options.discretization.Su_ux;
        Su_uy = options.discretization.Su_uy;
        Sv_vx = options.discretization.Sv_vx;
        Sv_vy = options.discretization.Sv_vy;
        
        ySu_ux = options.discretization.ySu_ux;
        ySu_uy = options.discretization.ySu_uy;
        ySv_vx = options.discretization.ySv_vx;
        ySv_vy = options.discretization.ySv_vy;
        
        dudx2 = (Su_ux*uh + ySu_ux).^2; % at the pressure points
        dudy2 = (Su_uy*uh + ySu_uy).^2;
        dvdx2 = (Sv_vx*vh + ySv_vx).^2;
        dvdy2 = (Sv_vy*vh + ySv_vy).^2;
        
        
        %% take weighted average to get contribution to a certain velocity
        % point
        weight = 1/2;
        
        % average dudx^2, effectively at u point
        % multiply by finite volume size
        diag1       = weight*ones(Nux_t,1);
        A1D         = spdiags([diag1 diag1],[0 1],Nux_t-2,Nux_t-1);
        Aux_u       = spdiags(Omu,0,Nu,Nu)*kron(speye(Nuy_in),A1D);
        dudx2_u     = Aux_u*dudx2;
        
        % average dudy^2, effectively at u point
        % multiply by finite volume size
        diag1       = weight*ones(Nuy_t,1);
        A1D         = spdiags([diag1 diag1],[0 1],Nuy_t-2,Nuy_t-1);
        Auy_u       = spdiags(Omu,0,Nu,Nu)*kron(A1D,speye(Nux_in));
        dudy2_u     = Auy_u*dudy2;
        
        % average dvdx^2, effectively at v point
        % multiply by finite volume size
        diag1       = weight*ones(Nvx_t,1);
        A1D         = spdiags([diag1 diag1],[0 1],Nvx_t-2,Nvx_t-1);
        Avx_v       = spdiags(Omv,0,Nv,Nv)*kron(speye(Nvy_in),A1D);
        dvdx2_v     = Avx_v*dvdx2;
        
        % average dvdy^2, effectively at v point
        % multiply by finite volume size
        diag1       = weight*ones(Nvy_t,1);
        A1D         = spdiags([diag1 diag1],[0 1],Nvy_t-2,Nvy_t-1);
        Avy_v       = spdiags(Omv,0,Nv,Nv)*kron(A1D,speye(Nvx_in));
        dvdy2_v     = Avy_v*dvdy2;
        
        
        % need to correct these dissipation terms because on "aligned" boundaries (d2udx2 and
        % d2vdy2) we do not get a proper summation by parts identity
        % "left boundary"
        switch BC.u.left
            case 'dir'
                dudx2_u(1:Nux_in:end) = dudx2_u(1:Nux_in:Nu) +  ...
                    0.5*hy.*(uh(1:Nux_in:Nu).^2)/hx(1);
        end
        
        switch BC.u.right
            case 'dir'
                % "right boundary"
                dudx2_u(Nux_in:Nux_in:end) = dudx2_u(Nux_in:Nux_in:Nu) + ...
                    0.5*hy.*(uh(Nux_in:Nux_in:Nu).^2)/hx(end);
        end
        
        switch BC.v.low
            case 'dir'
                % "bottom boundary"
                dvdy2_v(1:Nvx_in) = dvdy2_v(1:Nvx_in) + ...
                    0.5*hx.*(vh(1:Nvx_in).^2)/hy(1);
        end
        
        switch BC.v.up
            case 'dir'
                % "top boundary"
                dvdy2_v(Nv-Nvx_in+1:end) = dvdy2_v(Nv-Nvx_in+1:end) + ...
                    0.5*hx.*(vh(Nv-Nvx_in+1:end).^2)/hy(end);
        end
        
        %% then use definition of local kinetic energy to get Phi
        % since the local KE  is defined on a pressure volume, as sum of
        % its neighbouring velocity values (square), we get an additional
        % interpolation step
        
        % we need boundary contributions corresponding to ub^2, but we ignore these for now
        % average from velocity point to pressure point
        diag1     = weight*ones(Nux_t,1);
        A1D       = spdiags([diag1 diag1],[0 1],Npx,Npx+1);
        Au_k_BC  = ...
            BC_general(Npx+1,Nux_in,Npx+1-Nux_in,BC.u.left,BC.u.right,hx(1),hx(end));
        % extend to 2D
        Au_k      = kron(speye(Ny),A1D*Au_k_BC.B1D);
        
        % we need boundary contributions for vb^2, ignore these for now
        % average from velocity point to pressure point
        diag1     = weight*ones(Nvy_t,1);
        A1D       = spdiags([diag1 diag1],[0 1],Npy,Npy+1);
        Av_k_BC  = ...
            BC_general(Npy+1,Nvy_in,Npy+1-Nvy_in,BC.v.low,BC.v.up,hy(1),hy(end));
        % extend to 2D
        Av_k      = kron(A1D*Av_k_BC.B1D,speye(Nx));
        
        Phi = Au_k*(dudx2_u + dudy2_u) + Av_k*(dvdx2_v + dvdy2_v);
        
        
        % scale with the Reynolds number
        switch options.case.boussinesq
            
            case 'temp'
                % get T at v-locations
                Re = sqrt(options.temp.Ra/options.temp.Pr);
                
            otherwise
                Re = options.fluid.Re;
                
        end
        
        %
        Phi = (1/Re)*Phi;
        
        
        %% test correctness of dissipation
        % we compare integral of Phi to integral of V*(D*V)
        
        test = uh'*(options.discretization.Diffu*uh) + vh'*(options.discretization.Diffv*vh);
        
        if (abs(test + sum(sum(Phi)))>1e-14) % plus sign because Phi and V*D*V should have reverse signs
            warning('dissipation not consistent with V^T * DiffV * V; might be due to use of non-uniform grid, please check dissipation.m');
        end
        
        
        
        %% old tests used to construct the boundary contributions
        %         % note that Diffu is in integrated (finite volume) form, so the stencil is of the form
        %         % (dy/dx)*(1 -2 1)
        %         Diffu  = options.discretization.Dux*( options.discretization.Su_ux) + ...
        %             options.discretization.Duy*( options.discretization.Su_uy);
        %         Diffv  = options.discretization.Dvx*( options.discretization.Sv_vx) + ...
        %             options.discretization.Dvy*( options.discretization.Sv_vy);
        %
        %         testu = uh.*(Diffu*uh);
        %         testv = vh.*(Diffv*vh);
        %
        %         sum(testu)
        %         sum(testv)
        %
        %         sum(dudx2_u) + sum(dudy2_u)
        % %           0.5*sum(hy.*(uh(1:Nux_in:end).^2)/hx(1)) + ...
        % %           0.5*sum(hy.*(uh(Nux_in:Nux_in:end).^2)/hx(end))
        %         sum(dvdx2_v) + sum(dvdy2_v)
        % %           0.5*sum(hx.*(vh(1:Nvx_in).^2)/hy(1)) + ...
        % %           0.5*sum(hx.*(vh(end-Nvx_in+1:end).^2)/hy(end))
        
        
        
    otherwise
        
        error('dissipation function of specified viscosity model not implemented yet');
        
end


end

