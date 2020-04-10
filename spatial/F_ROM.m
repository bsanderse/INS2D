function [maxres,Fres,dF] = F_ROM(R,~,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
if (nargin<5)
    getJacobian = 0;
end

M  = options.rom.M;
B  = options.rom.B;

Om_inv = options.grid.Om_inv;

% steady inhomogeneous BC
Vbc = options.rom.Vbc;

% FOM velocity field (only needed when not precomputing)
if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0)
    V = B*R + Vbc;
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
    if (options.rom.weighted_norm == 0)
        conv  = B'*(Om_inv.*[convu;convv]);
        dconv = B'*(Om_inv.*[dconvu;dconvv])*B;
    elseif (options.rom.weighted_norm == 1)
        conv  = B'*[convu;convv];
        dconv = B'*[dconvu;dconvv]*B;
    end
    
end

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    [d2, dDiff] = diffusionROM(R,t,options,getJacobian);
elseif (options.rom.precompute_diffusion == 0)
    % approach 2: evaluate convection on FOM level, then map back
    [d2u,d2v,dDiffu,dDiffv] = diffusion(V,t,options,getJacobian);
    if (options.rom.weighted_norm == 0)   
        d2    = B'*(Om_inv.*[d2u;d2v]);
        dDiff = B'*(Om_inv.*[dDiffu;dDiffv])*B;
    elseif (options.rom.weighted_norm == 1)
        d2    = B'*[d2u;d2v];
        dDiff = B'*[dDiffu;dDiffv]*B;
    end
end

% body force
if (options.rom.precompute_force == 1)
    error('precomputed forcing term not implemented');
else
    [Fx, Fy] = force(t,options);
    if (options.rom.weighted_norm == 0)   
        F  = B'*(Om_inv.*[Fx;Fy]);
    elseif (options.rom.weighted_norm == 1)
        F  = B'*[Fx;Fy];
    end
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
    dF = spalloc(2*M,2*M,0);
end

end