function [maxres,Fres,dFres] = F_ROM(R,q,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
% inputs: R = ROM coefficients of velocity field; V = B*R
%         q = ROM coefficients of pressure field; p = Bp*q

if (nargin<5)
    getJacobian = 0;
end

if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0 || ...
    options.rom.precompute_force == 0)
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
    options.rom.precompute_force == 0)
    V = getFOM_velocity(R,t,options);
end

% FOM pressure field (only needed when not precomputing)
if (options.rom.div_free == 0 && options.rom.precompute_pressure == 0)
    p = getFOM_pressure(q,t,options);
end

% unsteady BC
if (options.BC.BC_unsteady == 1)
    if (options.rom.precompute_convection == 0 && options.rom.precompute_diffusion == 0)
        options = set_bc_vectors(t,options);
    else
        error('unsteady BC with precomputing not fully tested');
    end
end

% convection:
if (options.rom.precompute_convection == 1)
    % approach 1: (with precomputed matrices)
    [conv, dconv] = convectionROM(R,t,options,getJacobian);
elseif (options.rom.precompute_convection == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [convu, convv, dconvu, dconvv] = convection(V,V,t,options,getJacobian);
    conv  = B'*(Diag.*[convu;convv]);
    if (getJacobian == 1)
        dconv = B'*(Diag.*[dconvu;dconvv])*B;    
    end
end

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    [Diff, dDiff] = diffusionROM(R,t,options,getJacobian);
elseif (options.rom.precompute_diffusion == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [d2u,d2v,dDiffu,dDiffv] = diffusion(V,t,options,getJacobian);
    Diff  = B'*(Diag.*[d2u;d2v]);
    if (getJacobian == 1)
        dDiff = B'*(Diag.*[dDiffu;dDiffv])*B;
    end
end

% pressure gradient
% only needed if basis is not div-free
% note: this term is currently added in time_ERK_ROM directly
% if (options.rom.div_free == 0)
%     if (options.rom.precompute_pressure == 1)
%         % approach 1: (with precomputed matrices)
%         Gp = options.rom.G*q;
%     elseif (options.rom.precompute_pressure == 0)
%         % approach 2: evaluate convection on FOM level, then map back
%         Gp = options.rom.B' * options.discretization.G * p;
%     end
% end

% body force 
if (options.force.isforce == 1)
    if (options.rom.precompute_force == 1)
        F = options.rom.F;
        % this is a bit of a hack for the actuator ROM case with time dependent
        % body force, which prevents computing the projection of the force at
        % each time step
        
        if (options.force.force_unsteady == 1)               
            if (strcmp(options.case.project,'actuator_ROM'))
                if( options.time.t_start == t)
                    % give warning only once
                    warning('scaling unsteady force - actuator disk test case only!');
                end
                F = F*(1+sin(pi*t)); % see also pressure_additional_solve_ROM.m!
            else
                error('you have unsteady forcing with precomputation: check if settings are correct');
            end
        end
       
        if (getJacobian == 1)
            % Jacobian is not straightforward for general non-linear forcing    
            warning('precomputing Jacobian of force not available, using zero Jacobian');
            M  = options.rom.M;
            dF = spalloc(M,M,0);
        end
    else % no precomputing, use FOM expression and project to ROM (expensive)
        [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian);
        F   = B'*(Diag.*[Fx;Fy]);
        if (getJacobian == 1)
            dF  = B'*(Diag.*[dFx;dFy])*B;
        end
    %     dFy = spalloc(Nv,Nu+Nv,0);    
    end
else
    F  = options.rom.F;
    if (getJacobian == 1)
        M  = options.rom.M;
        dF = spalloc(M,M,0);
    end
end
% else
%     F = options.rom.F;
%     if (getJacobian == 1)
%         M  = options.rom.M;
%         dF = spalloc(M,M,0);
%     end
% end

% residual of ROM
Fres    = - conv + Diff + F;

% for the case of a non div-free basis, we have to add the pressure
% gradient
% note: this term is currently added in time_ERK_ROM directly
% if (options.rom.div_free == 0)
%     Fres = Fres - Gp;
% end

% for the case of non-orthogonal basis, we have to multiply with the
% inverse of (B'*B)
switch options.rom.rom_type
    
    case {'POD', 'Fourier'}
        
        
    case 'FDG' 
         
        % non-orthogonal basis, pre-multiply with (B^T*B)^{-1}
        Fres = options.rom.B_inv \ Fres;

    case 'FDG-Fourier'
        
        Fres = options.rom.B_inv .* Fres; 
        
    otherwise
        error ('wrong ROM type')
end

maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
       
    dFres   = -dconv + dDiff + dF;
    %     if (options.case.steady==0) % unsteady case, solve for velocities
    %         dF = spdiags(Om_inv,0,NV,NV)*dF;
    %     end
    
    switch options.rom.rom_type
        
        case {'POD', 'Fourier'}
            
            
        case 'FDG'
            
            dFres = options.rom.B_inv \ dFres;
            
        otherwise
            error ('wrong ROM type')
    end

else
    
    dFres = 0;
%     M  = options.rom.M;
%     dF = spalloc(2*M,2*M,0);
end

end