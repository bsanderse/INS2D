function [Phi, Jac] = dissipation(V,t,options,getJacobian)
% evaluate dissipation terms and optionally Jacobian
% the dissipation is defined here as (1/Re)* (grad u)^2
% the dissipation terms are defined on the pressure grid (where the kinetic
% energy is defined) and are a function of the velocity field V

% note the similarity with strain_tensor.m

visc = options.case.visc;
Np = options.grid.Np;
Nu = options.grid.Nu;

Nv = options.grid.Nv;

% Omu = options.grid.Omu;
% Omv = options.grid.Omv;


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
        
        dudx = (Su_ux*uh + ySu_ux);
        dudy = (Su_uy*uh + ySu_uy);
        dvdx = (Sv_vx*vh + ySv_vx);
        dvdy = (Sv_vy*vh + ySv_vy);
        
        %         dudx2 = (Su_ux*uh + ySu_ux).^2; % at the pressure points
        %         dudy2 = (Su_uy*uh + ySu_uy).^2;
        %         dvdx2 = (Sv_vx*vh + ySv_vx).^2;
        %         dvdy2 = (Sv_vy*vh + ySv_vy).^2;
        
        
        %         %% take weighted average to get contribution to a certain velocity
        %         % point
        %         weight = 1/2;
        %
        %         % average dudx^2, effectively at u point
        %         % multiply by finite volume size
        %         diag1       = weight*ones(Nux_t,1);
        %         A1D         = spdiags([diag1 diag1],[0 1],Nux_t-2,Nux_t-1);
        %
        %         % need to correct these dissipation terms because on "aligned" boundaries (d2udx2 and
        %         % d2vdy2) we do not get a proper summation by parts identity
        %         switch BC.u.left
        %             case 'dir'
        %                 % "left boundary"
        %                 A1D(1,1) = 1;
        % %                 dudx2_u(1:Nux_in:end) = dudx2_u(1:Nux_in:Nu) +  ...
        % %                     0.5*hy.*(uh(1:Nux_in:Nu).^2)/hx(1);
        %         end
        %         switch BC.u.right
        %             case 'dir'
        %                 % "right boundary"
        %                 A1D(end,end) = 1;
        % %                 dudx2_u(Nux_in:Nux_in:end) = dudx2_u(Nux_in:Nux_in:Nu) + ...
        % %                     0.5*hy.*(uh(Nux_in:Nux_in:Nu).^2)/hx(end);
        %         end
        %         Aux_u       = spdiags(Omu,0,Nu,Nu)*kron(speye(Nuy_in),A1D);
        %         dudx2_u     = Aux_u*(dudx.^2);
        %
        %         % average dudy^2, effectively at u point
        %         % multiply by finite volume size
        %         diag1       = weight*ones(Nuy_t,1);
        %         A1D         = spdiags([diag1 diag1],[0 1],Nuy_t-2,Nuy_t-1);
        %         Auy_u       = spdiags(Omu,0,Nu,Nu)*kron(A1D,speye(Nux_in));
        %         dudy2_u     = Auy_u*(dudy.^2);
        %
        %         % average dvdx^2, effectively at v point
        %         % multiply by finite volume size
        %         diag1       = weight*ones(Nvx_t,1);
        %         A1D         = spdiags([diag1 diag1],[0 1],Nvx_t-2,Nvx_t-1);
        %         Avx_v       = spdiags(Omv,0,Nv,Nv)*kron(speye(Nvy_in),A1D);
        %         dvdx2_v     = Avx_v*(dvdx.^2);
        %
        %         % average dvdy^2, effectively at v point
        %         % multiply by finite volume size
        %         diag1       = weight*ones(Nvy_t,1);
        %         A1D         = spdiags([diag1 diag1],[0 1],Nvy_t-2,Nvy_t-1);
        %
        %         switch BC.v.low
        %             case 'dir'
        %                 % "bottom boundary"
        %                 A1D(1,1) = 1;
        % %                 dvdy2_v(1:Nvx_in) = dvdy2_v(1:Nvx_in) + ...
        % %                     0.5*hx.*(vh(1:Nvx_in).^2)/hy(1);
        %         end
        %         switch BC.v.up
        %             case 'dir'
        %                 % "top boundary"
        %                 A1D(end,end) = 1;
        % %                 dvdy2_v(Nv-Nvx_in+1:end) = dvdy2_v(Nv-Nvx_in+1:end) + ...
        % %                     0.5*hx.*(vh(Nv-Nvx_in+1:end).^2)/hy(end);
        %         end
        %         Avy_v       = spdiags(Omv,0,Nv,Nv)*kron(A1D,speye(Nvx_in));
        %         dvdy2_v     = Avy_v*(dvdy.^2);
        %
        %
        %         %% then use definition of local kinetic energy to get Phi
        %         % since the local KE  is defined on a pressure volume, as sum of
        %         % its neighbouring velocity values (square), we get an additional
        %         % interpolation step
        %
        %         % we need boundary contributions corresponding to ub^2, but we ignore these for now
        %         % average from velocity point to pressure point
        %         diag1     = weight*ones(Nux_t,1);
        %         A1D       = spdiags([diag1 diag1],[0 1],Npx,Npx+1);
        %         Au_k_BC  = ...
        %             BC_general(Npx+1,Nux_in,Npx+1-Nux_in,BC.u.left,BC.u.right,hx(1),hx(end));
        %         % extend to 2D
        %         Au_k      = kron(speye(Ny),A1D*Au_k_BC.B1D);
        %
        %         % we need boundary contributions for vb^2, ignore these for now
        %         % average from velocity point to pressure point
        %         diag1     = weight*ones(Nvy_t,1);
        %         A1D       = spdiags([diag1 diag1],[0 1],Npy,Npy+1);
        %         Av_k_BC  = ...
        %             BC_general(Npy+1,Nvy_in,Npy+1-Nvy_in,BC.v.low,BC.v.up,hy(1),hy(end));
        %         % extend to 2D
        %         Av_k      = kron(A1D*Av_k_BC.B1D,speye(Nx));
        %
        
        Aux_u = options.discretization.Aux_u;
        Auy_u = options.discretization.Auy_u;
        Avx_v = options.discretization.Avx_v;
        Avy_v = options.discretization.Avy_v;
        Au_k = options.discretization.Au_k;
        Av_k = options.discretization.Av_k;
        
        dudx2_u     = Aux_u*(dudx.^2);
        dudy2_u     = Auy_u*(dudy.^2);
        dvdx2_v     = Avx_v*(dvdx.^2);
        dvdy2_v     = Avy_v*(dvdy.^2);
        
        Phi = Au_k*(dudx2_u + dudy2_u) + Av_k*(dvdx2_v + dvdy2_v);
        
        
        % scale with the Reynolds number or the 'effective' Reynolds number
        switch options.case.boussinesq
            
            case 'temp'
                %                 Re = sqrt(options.temp.Ra/options.temp.Pr);
                scale = options.temp.alfa1;
            otherwise
                scale = 1/options.fluid.Re;
        end
        
        Phi = scale*Phi;
        
        %
        
        
        %% test correctness of dissipation
        % we compare integral of Phi to integral of V*(D*V), where D
        % already includes the factor (1/Re)
        
        test = uh'*(options.discretization.Diffu*uh) + vh'*(options.discretization.Diffv*vh);
        
        if (abs(test + sum(Phi))/abs(test)>1e-12) % plus sign because Phi and V*D*V should have reverse signs
            warning('dissipation not consistent with V^T * DiffV * V; might be due to use of non-uniform grid or non-homogeneous BC; please check dissipation.m');
        end
        
        
        %% Jacobian
        if (getJacobian==1)
            
            N1 = length(dudx); %options.grid.N1;
            N2 = length(dudy); %options.grid.N2;
            N3 = length(dvdx); %options.grid.N3;
            N4 = length(dvdy); %options.grid.N4;
            
            Phi_u      = 2*Au_k*(Aux_u*spdiags(dudx,0,N1,N1)*Su_ux + Auy_u*spdiags(dudy,0,N2,N2)*Su_uy);
            Phi_v      = 2*Av_k*(Avx_v*spdiags(dvdx,0,N3,N3)*Sv_vx + Avy_v*spdiags(dvdy,0,N4,N4)*Sv_vy);
            
            Jac        = scale*[Phi_u Phi_v];
            
            
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

