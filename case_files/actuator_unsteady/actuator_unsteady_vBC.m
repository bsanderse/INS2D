function v = actuator_unsteady_vBC(x,y,t,options)
% boundary conditions for v for actuator

% coordinate left side domain:
x1 = options.grid.x1;

if ( (length(x)==1 && abs(x-x1)<1e-10)) % || (length(y)==1 && abs(y-0)<1e-10) )
    f     = 0.5;
    alpha = (pi/6)*sin(f*t);
    v     = sin(alpha)*ones(length(x)*length(y),1);
else    
    v = zeros(length(x)*length(y),1);    
end


end