%% mesh 
options = operator_mesh(options);


%% averaging operators
options = operator_averaging(options);


%% interpolation operators
options = operator_interpolation(options);


%% divergence (u,v)-> p and gradient p->(u,v) operator 
options = operator_divergence(options);


%% convection operators on u- and v- centered volumes
options = operator_convection_diffusion(options);


%% trbulence

% regularization modelling - this changes the convective term
if (regularize>0)
     
    operator_regularization;
    
end

% classical turbulence modelling via the diffusive term
switch visc

    case 'keps'   

        %% averaging viscosity to cell faces of u- and v- volumes
        ke_viscosity;

        %% k-e operators
        ke_convection;
        ke_diffusion;
        ke_production;

    case {'qr','LES','ML'}
        
        options = operator_turbulent_diffusion(options);
    
    case 'laminar'
        % do nothing, these are constructed in operator_convection_diffusion
        
    otherwise
        error('wrong value for visc parameter');
        
end


%% post-processing
options = operator_postprocessing(options);
