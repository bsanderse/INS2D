function u = actuator_Maarten_uBC(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;

if (length(x)==1 && abs(x-x1)<1e-10)
    u = options.fluid.U1*ones(length(y),1);
else
    u = zeros(length(x)*length(y),1);
end

end