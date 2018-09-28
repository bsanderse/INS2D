% constants for k-e model

% standard model
Ce1 = 1.44;
Ce2 = 1.92;
Cmu = 0.09;
sigmak  = 1.0;
sigmae  = 1.11;% 1.3; (1.3: standard value, 1.11: Richards and Hoxey)

% Crespo's model:
% Ce1 = 1.176;
% Ce2 = 1.92;
% Cmu = 0.033;
% sk  = 1.0;
% se  = 1.3;

% for MMS of Eca
% sigma = 4;
% sigma_nu = 2.5*sigma;
% kmax  = 0.01;
% numax = nu*1e3;


% inflow values for fully developed flow in NEUTRAL ABL

% density of air:
rho_inf = 1.225; % [kg/m^3]

% average roughness length at the mast3 location (from Jeroen)
% z0   = 0.017; % [m] 
% z0 = 2E-2; % rounded off
% z0 = 0.049; % @Sexbierum
z0 = 0.01; % Hargreaves

% Reference values (at the hub)
% z_ref = 80; % [m]
% U_ref = 10;    % [m/s]
% z_ref = 35; % [m] % @Sexbierum % also assume D = Hub height z_ref
% U_ref = 8;    % [m/s] % @Sexbierum

% Hargreaves:
z_ref = 6;
U_ref = 10;

z0_nondim = z0/z_ref;

% von Karman constant
kappa = 0.40;  % [-]

% friction velocity
u_fr  = kappa*U_ref/log((z0+z_ref)/z0);

% Model constant C_mu
C_mu  = 0.09; % Cmu from contants_ke.m

% Turbulent kinetic energy at ref. height (hub height)
TKE_inflow  = u_fr^2/sqrt(C_mu);    % [m^2/s^2]

% Turbulent (eddy) viscosity at ref. height (hub height)
epsilon_inf = u_fr^3/(kappa*z_ref); % [m^2/s^3]

%non-dimensionalization:
TKE_inflow  = TKE_inflow/U_ref^2;
epsilon_inf = epsilon_inf*z_ref/U_ref^3;