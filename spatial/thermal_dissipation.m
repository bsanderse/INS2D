function epsilonT = thermal_dissipation(T,t,options)
% evaluate thermal dissipation term
% the dissipation is defined here as (dT/dx)^2 + (dT/dy)^2
% note the similarity with dissipation.m

STx = options.discretization.STx;
STy = options.discretization.STy;
ySTx = options.discretization.ySTx;
ySTy = options.discretization.ySTy;

Npx = options.grid.Npx;
% Npy = options.grid.Npy;

%% squared derivatives
% this is in finite difference form (no weighting with volume sizes)
% note that alfa4 is not incorporated in these operators, as that would lead
% to alfa4^2
dTdx = (STx*T + ySTx); % ~1/Delta x
dTdy = (STy*T + ySTy); % ~1/Delta y
dTdx2 = dTdx.^2; % ~1/Delta x^2
dTdy2 = dTdy.^2; % ~1/Delta y^2

%% BC corrections; could maybe also be incorporated more systematically with volume weighting
% e.g. epsilonT = sum(options.grid.Omu*dTdx2) + sum(options.grid.Omv*dTdy2);
dTdx2_corr = dTdx2;
dTdy2_corr = dTdy2;

% dTdx left and right
dTdx2_corr(1:(Npx+1):end)       = 0.5*dTdx2(1:(Npx+1):end);
dTdx2_corr((Npx+1):(Npx+1):end) = 0.5*dTdx2((Npx+1):(Npx+1):end);
% dTdy up and low
dTdy2_corr(1:Npx)           = 0.5*dTdy2(1:Npx);
dTdy2_corr(end-Npx+1:end)   = 0.5*dTdy2(end-Npx+1:end);


%% scale by alfa4 and volume size
epsilonT = options.temp.alfa4*(sum(dTdx2_corr) + sum(dTdy2_corr))*options.grid.Omp(1); % dimensions dy/dx, dx/dy


%% test correctness of dissipation
% we compare epsilonT to T'*(Diff*T + yT);
% these should be the same but ONLY for homogeneous Neumann or homogeneous Dirichlet
% or periodic BC

% for nonhomogeneous BC, e.g. Dirichlet we get an additional term that
% will evaluates to something like the Nusselt number 

% general i: evalute T * Diff*T
% T_{i}*((T_{i+1} - T_{i}) - (T_{i} - T_{i-1})) = 
%  0.5*(T_{i+1} + T_{i})*(T_{i+1} - T_{i}) - 0.5*(T_{i+1} - T_{i})^2
% -0.5*(T_{i} + T_{i-1})*(T_{i} - T_{i-1}) - 0.5*(T_{i} - T_{i-1})^2

% i=1: T_{0} + T_{1} = 2*T_{b}, so T_{0} = 2*T_{b} - T_{1}, and T_{1} - T_{0} = 2*T_{1} - 2*ub
% T_{1}*((T_{2} - T_{1}) - (T_{1} - T_{0})) = 
%  0.5*(T_{2} + T_{1})*(T_{2} - T_{1}) - 0.5*(T_{2} - T_{1})^2
% -0.5*(T_{1} + T_{0})*(T_{1} - T_{0}) - 0.5*(T_{1} - T_{0})^2 =
%  0.5*(T_{2} + T_{1})*(T_{2} - T_{1}) - 0.5*(T_{2} - T_{1})^2
%               -2*T_{b}*(T_{1} - T_{b}) - 0.5*(2*T_{1} - 2*T_{b})^2

%  => for T_{b} = 0
%  0.5*(T_{2} + T_{1})*(T_{2} - T_{1}) - 0.5*(T_{2} - T_{1})^2
% -0.5*( (T_{1} - T_{b})/0.5)^2

%  => for (du/dx)_b = 0 (T_{1} - T_{0} = 0 or T_{b} = T_{1})
%  0.5*(T_{2} + T_{1})*(T_{2} - T_{1}) - 0.5*(T_{2} - T_{1})^2

% note that DiffT includes alfa4, but STx does not include alfa4
% note that DiffT includes volume sizes Omega, so we have dimensions dy/dx,
% dx/dy

test = T'*(options.discretization.DiffT*T + options.discretization.yDiffT);

if (abs(test + sum(epsilonT))>1e-14) % plus sign because epsilonT and T*D*T should have reverse signs
    warning('thermal dissipation not consistent with T^T * (DiffT * T + yT); might be due to use of non-uniform grid or non-homogeneous BC; please check thermal_dissipation.m');
end



end

