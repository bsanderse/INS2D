function options = constants_boussinesq(options)
% create the coefficients alfa_1...alfa_4 which are used in the Boussinesq equations
% the value of alfa depends on the choice of reference velocity in the
% non-dimensionalization

% mom. equation: Omega* dV/dt = -C(V) - G*p + alfa_1*DiffV + alfa2*T*e_y
% int. energy equation: Omega*dT/dt = -C_T(V,T) + alfa_3*Phi + alfa4*DiffT

Ra = options.temp.Ra;
Pr = options.temp.Pr;
Ge = options.temp.Ge;


switch options.temp.nondim_type
    
    case 1 % free fall velocity, uref = sqrt(beta*g*Delta T*H)
        alfa1 = sqrt(Pr/Ra);
        alfa2 = 1;
        alfa3 = Ge*sqrt(Pr/Ra);
        alfa4 = 1/sqrt(Pr*Ra);
        
    case 2 % uref = kappa/H (based on heat conduction time scale)
        alfa1 = Pr;
        alfa2 = Pr*Ra;
        alfa3 = Ge/Ra;
        alfa4 = 1;
        
    case 3 % uref = sqrt(c*Delta T)
        alfa1 = sqrt(Pr*Ge/Ra);
        alfa2 = Ge;
        alfa3 = sqrt(Pr*Ge/Ra);
        alfa4 = sqrt(Ge/(Pr*Ra));
        
        
    otherwise
        
        error('wrong choice for nondim_type in Boussinesq model');
        
end

% ratio between diffusive terms in momentum and dissipation in internal
% energy
gamma = alfa1/alfa3;

options.temp.alfa1 = alfa1;
options.temp.alfa2 = alfa2;
options.temp.alfa3 = alfa3;
options.temp.alfa4 = alfa4;
options.temp.gamma = gamma;