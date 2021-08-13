function [maxres,Fres,dFres] = F_ROM(R,~,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field

% warning pressure contribution for outflow BC not included
if (nargin<5)
    getJacobian = 0;
end

%% pfusch
% if (options.rom.precompute_convection == 0 && options.rom.precompute_diffusion == 0 && ...
%     options.rom.precompute_force == 0)
%     V = getFOM_velocity(R,t,options);
%     B  = options.rom.B;
% %     Diag = options.grid.Om_inv; eba net
%     nopressure = 1;
%     [maxres,Fres,dFres] = myF(V,V,0,t,options,getJacobian,nopressure);
% %     Fres2 = B'*(Diag.*Fres); eba net
%     Fres = B'*Fres;
%     maxres  = max(abs(Fres));
%     return
% end
%%

% if (options.rom.bc_recon == 2)
%     if (getJacobian == 1)
%         error('bc recon = 2 and get Jacobian = 1 not implemented')
%     end
%     V = getFOM_velocity(R,t,options);
%     p = zeros(options.grid.Nx*options.grid.Ny,1);
%     p = pressure_additional_solve(V,p,t,options);
%     nopressure = 0;
% %     cd spatial
% %     ls
%     [maxres,Fres,dFres] = myF(V,V,p,t,options,getJacobian,nopressure);
%     B = options.rom.B;
%     Diag = options.grid.Om_inv;
%     Fres = B'*(Diag.*Fres);
% %     cd ..
%     return
% end

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

% unsteady BC
if (options.BC.BC_unsteady == 1)
    if (options.rom.precompute_convection == 0 || options.rom.precompute_diffusion == 0)
        options = set_bc_vectors(t,options);
    elseif (options.rom.bc_recon ~= 1 && options.rom.bc_recon ~= 3)
        error('unsteady BC with precomputing not fully tested');
    end
end
%         options = set_bc_vectors(t,options);
% disp('F ROM 34 manipulated')

% convection:
if (options.rom.precompute_convection == 1)
    % approach 1: (with precomputed matrices)
    if (options.rom.rom_bc == 2)
        if options.rom.bc_recon == 3
%             [conv1, dconv1] = convectionROM_unsteadyBC2(R,t,options,getJacobian);
            [conv, dconv] = convectionROM_unsteadyBC2(R,t,options,getJacobian);
        else % debugging
            [conv, dconv] = convectionROM_unsteadyBC(R,t,options,getJacobian);
%         [conv0, dconv0] = convectionROM_unsteadyBC(R,t,options,getJacobian);
%         norm(conv0-conv)
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

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    if (options.rom.rom_bc == 2)
        if options.rom.bc_recon == 3
%             [Diff1, dDiff1] = diffusionROM_unsteadyBC2(R,t,options,getJacobian);
            [Diff, dDiff] = diffusionROM_unsteadyBC2(R,t,options,getJacobian);
        else % debugging
            [Diff, dDiff] = diffusionROM_unsteadyBC(R,t,options,getJacobian);
%         [Diff0, dDiff0] = diffusionROM_unsteadyBC(R,t,options,getJacobian);
%         norm(Diff0-Diff)
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
%     [d2u,d2v,dDiffu,dDiffv] = mydiffusion(V,t,options,getJacobian);
%     Diff1  = B'*(Diag.*[d2u;d2v]);
%     norm(Diff-Diff1)
%     Diff = Diff1;
% %     dDiff = B'*(Diag.*[dDiffu;dDiffv])*B;
% %     disp('F ROM 64 manipulated')


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
        warning('precomputing Jacobian of force not available, using zero Jacobian');
        M  = options.rom.M;
        dF = spalloc(M,M,0);
    end
else
    [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian);
    F   = B'*(Diag.*[Fx;Fy]);
    dF  = B'*(Diag.*[dFx;dFy])*B;
%     dFy = spalloc(Nv,Nu+Nv,0);    
end

% residual of ROM
Fres    = - conv + Diff + F;
% [d2u,d2v] = mydiffusion(get_unsteadyVbc(t,options),0,options,0);
% DiffBC = B'*[d2u;d2v];
% Fres    = - conv + Diff + F + DiffBC;

% norm(Fres-Fres2)

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
    
else
    
    dFres = 0;
%     M  = options.rom.M;
%     dF = spalloc(2*M,2*M,0);
end

end