% pressure correction steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Having set up the matrices Ru and Rv, the right-hand             %
%  side vector for the Poisson equation for the pressure is set up. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary condition for the difference in pressure between time
% steps; only non-zero in case of fluctuating outlet pressure
y_dp = zeros(Nu+Nv,1);

% divergence of Ru and Rv is directly calculated with M
f    = (1/dt)*(M*R + yM) - M*y_dp;

dp   = pressure_poisson(f,t,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Update the velocity field                                        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V   = R - dt*Om_inv.*(G*dp + y_dp);      