function [maxres,Fres,dF] = FTemp(T,V,t,options,getJacobian)
% calculate rhs of temperature equations and, optionally, Jacobian with respect to velocity
% field and temperature
% T: temperature
% V: velocity field

% getJacobian = 1: return dFdT and dFdV
if (nargin<5)
    getJacobian = 0;
end

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = options.grid.NV;

% unsteady BC
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end


% temperature
switch options.case.boussinesq
    
    case 'temp'
        Fv   = Fv - T;
        FT   = conv_diff_temperature(T,V,t,options,getJacobian);
        
        Fres = [Fu;Fv;FT];
        
    otherwise
        Fres = [Fu;Fv];    
end

% norm of residual
maxres  = max(abs(Fres));

if (getJacobian==1)
    % Jacobian requested
    % we return only the Jacobian with respect to V (not p)
           
    dFu  = - dconvu + dDiffu + dFx;
    dFv  = - dconvv + dDiffv + dFy;
    
    dF   = [dFu; dFv];
    
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end