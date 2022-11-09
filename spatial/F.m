function [maxres,F,dF] = F(V,C,p,T,t,options,getJacobian,nopressure)
%,FTemp,dFTemp]
% calculate rhs of momentum equations and, optionally, Jacobian with respect to velocity
% field
% in addition, rhs of temperature equation is computed (FTemp), and its
% Jacobian with respect to temperature and velocity
% V: velocity field
% C: 'convection' field: e.g. d(c_x u)/dx + d(c_y u)/dy; usually c_x = u,
% c_y=v, so C=V
% p: pressure

% getJacobian = 1: return dFdV
% nopressure = 1: exclude pressure gradient; in this case input argument p is not used

if (nargin<8)
    nopressure = 0;    
end
if (nargin<7)
    getJacobian = 0;
end

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = options.grid.NV;

% unsteady BC
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end

if (nopressure == 0)
    % pressure
    Gx   = options.discretization.Gx;
    Gy   = options.discretization.Gy;
    y_px = options.discretization.y_px;
    y_py = options.discretization.y_py;
    Gpx  = Gx*p + y_px;
    Gpy  = Gy*p + y_py;
end

% convection:
[convu, convv, dconvu, dconvv] = convection(V,C,t,options,getJacobian);

% diffusion
[d2u, d2v, dDiffu, dDiffv] = diffusion(V,t,options,getJacobian);

% body force
if (options.force.isforce == 1)
    if (options.force.force_unsteady == 1)
        [Force_x, Force_y, dForce_x, dForce_y] = force(V,t,options,getJacobian);
    else
        Force_x = options.force.Fx;
        Force_y = options.force.Fy;
        dForce_x = spalloc(Nu,NV,0);
        dForce_y = spalloc(Nv,NV,0);            
    end
else
    Force_x  = zeros(Nu,1);
    Force_y  = zeros(Nv,1);
    dForce_x = spalloc(Nu,NV,0);
    dForce_y = spalloc(Nv,NV,0);      
end



% residual in Finite Volume form
Fu   = - convu + d2u + Force_x;
Fv   = - convv + d2v + Force_y;

% nopressure=0 is the most common situation, in which we return the entire 
% right-hand side vector
if (nopressure == 0) 
    %residual including the pressure contribution
    Fu = Fu - Gpx;
    Fv = Fv - Gpy;
end


% addition of temperature term in momentum equation 
switch options.case.boussinesq
    
    case 'temp'
        % get T at v-locations
        Fv     = Fv + options.discretization.AT_v*T; 
end

FV = [Fu;Fv];


% additional temperature equation
switch options.case.boussinesq
    
    case 'temp'
        [FTemp,dFTemp,dFTemp_V]  = conv_diff_temperature(T,V,t,options,getJacobian);

        switch options.temp.dissipation
            case 1
                %  add dissipation to internal energy equation
                Phi = dissipation(V,t,options,getJacobian);

                FTemp = FTemp - Phi;
        end
        
        F = [FV; FTemp];

    otherwise
        F = FV;
end

% norm of residual
maxres  = max(abs(F));

if (getJacobian==1)
    % Jacobian requested
    % we return only the Jacobian with respect to V (not p)
           
    dFu  = - dconvu + dDiffu + dForce_x;
    dFv  = - dconvv + dDiffv + dForce_y;
    

    % additional temperature equation
    switch options.case.boussinesq

        case 'temp'
            % add effect of temperature in v-momentum equations
            % add temperature equation
            NT = options.grid.NT;
            dF   = [dFu spalloc(Nu,NT,0); ...
                    dFv options.discretization.AT_v; ...
                    dFTemp_V dFTemp];
    
        otherwise
            dF   = [dFu; dFv];
    end
    
else
    dF = spalloc(Nu+Nv,Nu+Nv,0);
end

end