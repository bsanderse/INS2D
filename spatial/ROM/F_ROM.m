function [maxres,Fres,dF] = F_ROM(R,~,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
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
if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0)
    V = getFOM_velocity(R,t,options);
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
    dconv = B'*(Diag.*[dconvu;dconvv])*B;    
end

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    [d2, dDiff] = diffusionROM(R,t,options,getJacobian);
elseif (options.rom.precompute_diffusion == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [d2u,d2v,dDiffu,dDiffv] = diffusion(V,t,options,getJacobian);
    d2    = B'*(Diag.*[d2u;d2v]);
    dDiff = B'*(Diag.*[dDiffu;dDiffv])*B;
end

% body force
if (options.rom.precompute_force == 1)
    F = options.rom.F;
else
    [Fx, Fy] = force(t,options);
    F  = B'*(Diag.*[Fx;Fy]);
end

% residual of ROM
Fres    = - conv + d2 + F;


maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
       
    dF   = -dconv + dDiff;
    %     if (options.case.steady==0) % unsteady case, solve for velocities
    %         dF = spdiags(Om_inv,0,NV,NV)*dF;
    %     end
    
else
    M  = options.rom.M;
    dF = spalloc(2*M,2*M,0);
end

end