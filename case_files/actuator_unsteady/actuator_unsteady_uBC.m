function u = actuator_unsteady_uBC(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;
% coordinate lower side domain:
% y1 = options.grid.y1;


if ( (length(x)==1 && abs(x-x1)<1e-10)) %|| (length(y)==1 && abs(y-y1)<1e-10) )
    f     = 0.5;
    alpha = (pi/6)*sin(f*t);
    u     = cos(alpha)*ones(length(x)*length(y),1);
else
    u = zeros(length(x)*length(y),1);
end

end