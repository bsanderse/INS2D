function u = actuator_uBC(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;

if (length(x)==1 && abs(x-x1)<1e-10)
%     Ct = 0.001;
%     u = 1 + 0.5*Ct*1 / (2*pi) * (atan((1/2-y)/x) + atan((1/2+y)/x)) ;
    u = ones(length(y),1);
else
    u = zeros(length(x)*length(y),1);
end

end