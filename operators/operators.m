%% mesh 
operator_mesh;


%% averaging operators
operator_averaging;

%% interpolation operators
operator_interpolation;

%% divergence (u,v)-> p and gradient p->(u,v) operator 
operator_divergence;

%% convection and diffusion operators on u- and v- centered volumes
operator_convection_diffusion;





%% post-processing
operator_postprocessing;


%% boundary conditions
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;


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

