function [maxres,Fres,dFres] = F_ROM(R,~,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field

% warning pressure contribution for outflow BC not included
if (nargin<5)
    getJacobian = 0;
end

% check whether open boundary conditions occur
BC = options.BC;
obc = max(strcmp({BC.u.left,BC.u.right,BC.u.low,BC.u.up},'mvp-obc')); 


if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0 || ...
    options.rom.precompute_force == 0      || (obc && options.rom.precompute_obc == 0))
    B  = options.rom.B;
    if (options.rom.weighted_norm == 0)
        Diag = options.grid.Om_inv;
    elseif (options.rom.weighted_norm == 1)
        NV   = options.grid.Nu + options.grid.Nv;
        Diag = ones(NV,1);
    end
end
    
% FOM velocity field (only needed when not precomputing)
if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0 || ...
    options.rom.precompute_force == 0      || (obc && options.rom.precompute_obc == 0))
    V = getFOM_velocity(R,t,options);
end

% unsteady BC
if (options.BC.BC_unsteady == 1)
    if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0)
        options = set_bc_vectors(t,options);
    elseif (options.rom.bc_recon ~= 1 && options.rom.bc_recon ~= 3)
        error('unsteady BC with precomputing not fully tested');
    end
end

% convection:
if (options.rom.precompute_convection == 1)
    % approach 1: (with precomputed matrices)
    if (options.rom.rom_bc == 2)
        if options.rom.bc_recon == 3
            [conv, dconv] = convectionROM_unsteadyBC2(R,t,options,getJacobian);
        else 
            [conv, dconv] = convectionROM_unsteadyBC(R,t,options,getJacobian);
        end
    else
        [conv, dconv] = convectionROM(R,t,options,getJacobian);
    end
elseif (options.rom.precompute_convection == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [convu, convv, dconvu, dconvv] = convection(V,V,t,options,getJacobian);
    conv  = B'*(Diag.*[convu;convv]);
    dconv = B'*(Diag.*[dconvu;dconvv])*B;  
end
%% testing
%     options = set_bc_vectors(t,options);
%     V = getFOM_velocity(R,t,options);
%     B  = options.rom.B;
%     NV   = options.grid.Nu + options.grid.Nv;
%     Diag = ones(NV,1);
% 
%     [convu, convv, dconvu, dconvv] = convection(V,V,t,options,getJacobian);
%     conv_test  = B'*(Diag.*[convu;convv]);
%     norm(conv-conv_test)
%%

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    if (options.rom.rom_bc == 2)
        if options.rom.bc_recon == 3
            [Diff, dDiff] = diffusionROM_unsteadyBC2(R,t,options,getJacobian);
        else 
            [Diff, dDiff] = diffusionROM_unsteadyBC(R,t,options,getJacobian);
        end
    else
        [Diff, dDiff] = diffusionROM(R,t,options,getJacobian);
    end
elseif (options.rom.precompute_diffusion == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [d2u,d2v,dDiffu,dDiffv] = mydiffusion(V,t,options,getJacobian);
    Diff  = B'*(Diag.*[d2u;d2v]);
    dDiff = B'*(Diag.*[dDiffu;dDiffv])*B;
end
%% testing
%     [d2u,d2v,dDiffu,dDiffv] = mydiffusion(V,t,options,getJacobian);
%     Diff_test  = B'*(Diag.*[d2u;d2v]);
%     norm(Diff-Diff_test)
%%

% body force
if (options.rom.precompute_force == 1)
    F = options.rom.F;
    % this is a bit of a hack for the actuator ROM case with time dependent
    % body force, which prevents computing the projection of the force at
    % each time step
    if (options.case.force_unsteady == 1)
        F = F*(1+sin(pi*t)); % see also pressure_additional_solve_ROM.m!
    end
    if (getJacobian == 1)
        % Jacobian is not straightforward for general non-linear forcing 
%         botch:
%         warning('precomputing Jacobian of force not available, using zero Jacobian');
        M  = options.rom.M;
        dF = spalloc(M,M,0);
    end
else
    [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian);
    F   = B'*(Diag.*[Fx;Fy]);
    dF  = B'*(Diag.*[dFx;dFy]);
%     dFy = spalloc(Nv,Nu+Nv,0);    
end

%% testing
%     [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian);
%     F_test   = B'*(Diag.*[Fx;Fy]); 
%     norm(F-F_test)
%%

% open boundaries
if obc
    if (options.rom.precompute_obc == 1)
%         if (options.rom.rom_bc == 2)
%             if options.rom.bc_recon == 3
                [y_O_diag_ROM, Jac_y_O_diag_ROM] = obcROM(R,t,options,getJacobian);
                switch BC.gO_type
                    case 0
                        gO_ROM = 0;
                        Jac_gO_ROM = 0;
                    case 1
                        [gO_ROM, Jac_gO_ROM] = gOROM(R,t,options,getJacobian);
                    case 2
                        error('Sorry, obc offline decomposition for more complex gO not implemented')
                end
                y_O_ROM = y_O_diag_ROM - gO_ROM;
                if (getJacobian==1)
                    Jac_y_O_ROM = Jac_y_O_diag_ROM - Jac_gO_ROM;
                end
%             else
%                 error('Sorry, precomputation of obc not implemented for bc_recon=/=3')
%             end
%         else
%             error('Sorry, precomputation of obc not implemented for rom_bc=/=2')
%         end
    else
        gO = @(V) options.BC.gO(V) + 0*V;
        dgO = @(V) options.BC.dgO(V) + 0*V;
        
        Conv_diag = options.grid.C;
        y_O_diag = Conv_diag*V;
        
        id_normal = options.grid.id_normal;
        id_tangential = options.grid.id_tangential;
        id_n_t = id_normal+id_tangential;
        V_n_t = id_n_t.*V;
        
%         y_O2  = y_O_diag.*V - V_n_t.*gO(V);
        gO_factor = options.grid.gO_factor;
        y_O = spdiags(V,0,NV,NV)*(Conv_diag*V)...
            - spdiags(gO_factor,0,NV,NV)*(spdiags(gO(V),0,NV,NV)*V);
%         norm(y_O-y_O2)
        y_O_ROM = B'*(Diag.*y_O);
        if (getJacobian==1)
            gO_factor = options.grid.gO_factor;
            Jac_y_O = spdiags(Conv_diag*V,0,NV,NV) ...
            + spdiags(V,0,NV,NV)*Conv_diag ...
            - spdiags(gO_factor,0,NV,NV)*spdiags(gO(V),0,NV,NV) ...
            - spdiags(gO_factor,0,NV,NV)*spdiags(dgO(V),0,NV,NV) ...
              * spdiags(V,0,NV,NV);
          Jac_y_O_ROM = B'*(Diag.*Jac_y_O)*B; % just a guess
        end
    end
    %% testing
        V = getFOM_velocity(R,t,options);
        NV = options.grid.NV;
        B = options.rom.B;
            if (options.rom.weighted_norm == 0)
        Diag = options.grid.Om_inv;
    elseif (options.rom.weighted_norm == 1)
        NV   = options.grid.Nu + options.grid.Nv;
        Diag = ones(NV,1);
    end

            gO = @(V) options.BC.gO(V) + 0*V;
        dgO = @(V) options.BC.dgO(V) + 0*V;
        
        Conv_diag = options.grid.C;
        y_O_diag = Conv_diag*V;
        
        id_normal = options.grid.id_normal;
        id_tangential = options.grid.id_tangential;
        id_n_t = id_normal+id_tangential;
        V_n_t = id_n_t.*V;
        
%         y_O2  = y_O_diag.*V - V_n_t.*gO(V);
        gO_factor = options.grid.gO_factor;
        y_O = spdiags(V,0,NV,NV)*(Conv_diag*V)...
            - spdiags(gO_factor,0,NV,NV)*(spdiags(gO(V),0,NV,NV)*V);
        y_O_ROM2 = B'*(Diag.*y_O);
        norm(y_O_ROM2-y_O_ROM)
        17
            
%         gO = options.BC.gO;
%         dgO = options.BC.dgO;
%         
%         Conv_diag = options.grid.C;
%         y_O_diag = Conv_diag*V;
%         
%         id_normal = options.grid.id_normal;
%         id_tangential = options.grid.id_tangential;
%         id_n_t = id_normal+id_tangential;
%         V_n_t = id_n_t.*V;
%         
%         y_O  = y_O_diag.*V - V_n_t.*gO(V);
%         y_O_ROM_test = B'*(Diag.*y_O);
%         norm(y_O_ROM-y_O_ROM_test)
    %%
end

% residual of ROM
Fres    = - conv + Diff + F;
if obc
    Fres = Fres + y_O_ROM;
end

maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
       
    dFres   = -dconv + dDiff + dF;
    if obc
        dFres = dFres + Jac_y_O_ROM;
    end
    %     if (options.case.steady==0) % unsteady case, solve for velocities
    %         dF = spdiags(Om_inv,0,NV,NV)*dF;
    %     end
    
else
    
    dFres = 0;
%     M  = options.rom.M;
%     dF = spalloc(2*M,2*M,0);
end

end