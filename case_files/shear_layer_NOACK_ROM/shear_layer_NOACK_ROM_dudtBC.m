function dudt = shear_layer_NOACK_ROM_dudtBC(x,y,t,options)
% boundary conditions for u for shear layer

% coordinate left side domain:
x1 = options.grid.x1;


if (length(x)==1 && abs(x-x1)<1e-10)

    
    mode_u = options.fluid.mode_u;
    omega  = options.fluid.omega;
    alpha  = options.fluid.alpha;    
    pert   = options.fluid.pert;
    
    I      = sqrt(-1);
    du_pert = real((-I*omega)*mode_u(y).*exp(I*(alpha*x - omega*t)));

    dudt    = pert*du_pert;
    
    
else
    dudt = zeros(length(x)*length(y),1);
end

end