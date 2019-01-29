function dvdt = actuator_unsteady_dvdtBC(x,y,t,options)
% boundary conditions for v for actuator

% coordinate left side domain:
x1 = options.grid.x1;

if ( (length(x)==1 && abs(x-x1)<1e-10)) % || (length(y)==1 && abs(y-0)<1e-10) )
    f     = 0.5;
    a     = (pi/6);
    alpha = a*sin(f*t);
    dvdt  = cos(alpha)*a*cos(f*t)*f*ones(length(x)*length(y),1);
else    
    dvdt  = zeros(length(x)*length(y),1);    
end


end