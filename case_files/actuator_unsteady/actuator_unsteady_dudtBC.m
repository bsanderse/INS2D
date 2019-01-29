function dudt = actuator_unsteady_dudtBC(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;
% coordinate lower side domain:
% y1 = options.grid.y1;


if ( (length(x)==1 && abs(x-x1)<1e-10)) %|| (length(y)==1 && abs(y-y1)<1e-10) )
    f     = 0.5;
    a     = (pi/6);
    alpha = a*sin(f*t);
%     u     = cos(alpha)*ones(length(x)*length(y),1);
%     a     = (pi/6);
%     f     = 0.5;
    dudt  = -sin(alpha)*a*cos(f*t)*f*ones(length(x)*length(y),1);    
else
    dudt = zeros(length(x)*length(y),1);
end

end