function v = shear_layer_NOACK_ROM_vBC(x,y,t,options)
% boundary conditions for v for shear layer

% coordinate left side domain:
x1 = options.grid.x1;
    
    if (length(x)==1 && abs(x-x1)<1e-10)

        mode_v = options.fluid.mode_v;
        omega  = options.fluid.omega;
        alpha  = options.fluid.alpha;
        pert   = options.fluid.pert;

        I      = sqrt(-1);

        v_pert = real(mode_v(y).*exp(I*(alpha*x - omega*t)));

        v = pert*v_pert;


    else
        v = zeros(length(x)*length(y),1);
    end
end