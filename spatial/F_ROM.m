function [maxres,Fres,dF] = F_ROM(R,p,t,options,getJacobian)
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

% ru = R(1:M);
% rv = R(M+1:2*M);

% unsteady BC
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end

% Gx   = options.discretization.Gx;
% Gy   = options.discretization.Gy;
% y_px = options.discretization.y_px;
% y_py = options.discretization.y_py;
Gpx  = zeros(M,1); %Gx*p + y_px;
Gpy  = zeros(M,1); %Gy*p + y_py;

% convection:
% approach 1: (with precomputed matrices) 
% [convu, convv, dconvu, dconvv] = convectionROM(ru,rv,t,options,getJacobian);
% approach 2:
[convu,convv] = convection(Bu*R,Bv*R,t,options,0);
% convu = Bu'*(options.grid.Omu_inv.*convu);
% convv = Bv'*(options.grid.Omv_inv.*convv);
% conv  = convu+convv;
conv  = B'*(options.grid.Om_inv.*[convu;convv]);

% diffusion
% approach 1: (with precomputed matrices) 
[d2, dDiff] = diffusionROM(R,t,options,0);
% approach 2:
% [d2u,d2v] = diffusion(Bu*R,Bv*R,t,options,0);
% d2 = B'*Om_inv.*[d2u;d2v];

% body force
% [Fx, Fy] = force(t,options);
% Fx = Bu'*(options.grid.Omu_inv.*Fx);
% Fy = Bv'*(options.grid.Omv_inv.*Fy);
% F  = Fx + Fy;

% residual of ROM
Fres    = - conv + d2;
% Fu   = - convu + d2u - Gpx + Fx;
% Fv   = - convv + d2v - Gpy + Fy;

% if (options.case.steady==0) % unsteady case, solve for velocities
%     Fres = Om_inv.*[Fu;Fv];
% else
%     Fres = [Fu;Fv];
% end

maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
        
    % we assume here that the body force Fx, Fy is not depending on the
    % solution u,v
    % so we only have convection and diffusion in the Jacobian
    
    dFu  = - dconvu + dDiffu;
    dFv  = - dconvv + dDiffv;
    
    dF   = [dFu; dFv];
    
%     if (options.case.steady==0) % unsteady case, solve for velocities
%         dF = spdiags(Om_inv,0,NV,NV)*dF;
%     end
    
else
    dF = spalloc(2*M,2*M,0);
end

end