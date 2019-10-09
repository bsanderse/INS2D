function u = shear_layer_NOACK_ROM_uBC(x,y,t,options)
% boundary conditions for u for shear layer

delta = options.fluid.d_layer;
U1    = options.fluid.U1;
U2    = options.fluid.U2;

% coordinate left side domain:
x1 = options.grid.x1;

y1 = options.grid.y1;
y2 = options.grid.y2;


if (length(x)==1 && abs(x-x1)<1e-10)

    % base flow + pert at inflow
    u = 0.5*(U1 + U2) + 0.5*(U1-U2)*tanh(y./delta);

    
    mode_u = options.fluid.mode_u;
    omega  = options.fluid.omega;
    alpha  = options.fluid.alpha;    
    pert   = options.fluid.pert;
    
    I      = sqrt(-1);
    u_pert = real(mode_u(y).*exp(I*(alpha*x - omega*t)));

    u = u + pert*u_pert;
    
elseif (length(y)==1 && abs(y-y1)<1e-10 && strcmp(options.BC.u.low,'dir'))

    % base flow at top and bottom
    u = (0.5*(U1 + U2) + 0.5*(U1-U2)*tanh(y./delta))*ones(length(x),1);
    
    
elseif (length(y)==1 && abs(y-y2)<1e-10 && strcmp(options.BC.u.up,'dir')) 
    
    % base flow at top and bottom
    u = (0.5*(U1 + U2) + 0.5*(U1-U2)*tanh(y./delta))*ones(length(x),1);
        
else
    u = zeros(length(x)*length(y),1);
end

end