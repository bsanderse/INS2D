%% mesh 
options = operator_mesh(options);


%% averaging operators
options = operator_averaging(options);


%% interpolation operators
options = operator_interpolation(options);


%% divergence (u,v)-> p and gradient p->(u,v) operator 
options = operator_divergence(options);


%% convection and diffusion operators on u- and v- centered volumes
options = operator_convection_diffusion(options);


%% post-processing
options = operator_postprocessing(options);


%% turbulence
if (regularize>0)
     
    operator_regularization;
    
end

if ( strcmp(visc,'turbulent') )

    %% averaging viscosity to cell faces of u- and v- volumes
    operators_viscosity;

    %% k-e operators
    ke_convection;
    ke_diffusion;
    ke_production;

end
