function u = uBC_(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;
% coordinate lower side domain:
% y1 = options.grid.y1;


if ( (length(x)==1 && abs(x-x1)<1e-10))
    u = ones(length(x)*length(y),1);
else
    u = zeros(length(x)*length(y),1);
end

end