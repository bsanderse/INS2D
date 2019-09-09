function [maxres,Fres,dF] = F_ROM(R,~,t,options,getJacobian)
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
if (nargin<5)
    getJacobian = 0;
end

M  = options.rom.M;
B  = options.rom.B;
Bu = options.rom.Bu;
Bv = options.rom.Bv;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;

Om_inv = options.grid.Om_inv;

% unsteady BC
if (options.BC.BC_unsteady == 1)
    if (options.rom.precompute_convection == 0 && options.rom.precompute_diffusion == 0)
        options = set_bc_vectors(t,options);
    else
        error('unsteady BC with precomputing not fully tested');
    end
end

% Gx   = options.discretization.Gx;
% Gy   = options.discretization.Gy;
% y_px = options.discretization.y_px;
% y_py = options.discretization.y_py;
Gpx  = zeros(M,1); %Gx*p + y_px;
Gpy  = zeros(M,1); %Gy*p + y_py;

% convection:
if (options.rom.precompute_convection == 1)
    % approach 1: (with precomputed matrices)
    error('precomputed convection term not fully tested');
%     [convu, convv, dconvu, dconvv] = convectionROM(ru,rv,t,options,getJacobian);
elseif (options.rom.precompute_convection == 0)
    % approach 2:
    [convu, convv, dconvu, dconvv] = convection(B*R,t,options,getJacobian);
    conv  = B'*(Om_inv.*[convu;convv]);
    dconv = B'*(Om_inv.*[dconvu;dconvv])*B;

end

% diffusion
if (options.rom.precompute_diffusion == 1)
    % approach 1: (with precomputed matrices)
    [d2, dDiff] = diffusionROM(R,t,options,getJacobian);
elseif (options.rom.precompute_diffusion == 0)
    % approach 2:
    [d2u,d2v,dDiffu,dDiffv] = diffusion(B*R,t,options,getJacobian);
    d2    = B'*(Om_inv.*[d2u;d2v]);
    dDiff = B'*(Om_inv.*[dDiffu;dDiffv]);
end

% body force
if (options.rom.precompute_force == 1)
    error('precomputed forcing term not implemented');
else
    [Fx, Fy] = force(t,options);
%     Fx = Bu'*(options.grid.Omu_inv.*Fx);
%     Fy = Bv'*(options.grid.Omv_inv.*Fy);
%     F  = Fx + Fy;
    F  = B'*(Om_inv.*[Fx;Fy]);
end

% residual of ROM
Fres    = - conv + d2 + F;
% Fu   = - convu + d2u - Gpx + Fx;
% Fv   = - convv + d2v - Gpy + Fy;


maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
        
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
    
%     dFu  = - dconvu + dDiffu;
%     dFv  = - dconvv + dDiffv;
%     
%     dF   = [dFu; dFv];
    
    dF   = -dconv + dDiff;
%     if (options.case.steady==0) % unsteady case, solve for velocities
%         dF = spdiags(Om_inv,0,NV,NV)*dF;
%     end
    
else
    dF = spalloc(2*M,2*M,0);
end

end