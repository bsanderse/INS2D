function u = mixing_layer_uBC(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;
% coordinate lower side domain:
% y1 = options.grid.y1;


if ( (length(x)==1 && abs(x-x1)<1e-10)) %|| (length(y)==1 && abs(y-y1)<1e-10) )
    Delta_U = 1;
    U_avg = 1;
    eps1 = 0.082*U_avg;
    eps2 = 0.018*U_avg;
    n1 = 0.4*pi;
    n2 = 0.3*pi;
    om1 = 0.22;
    om2 = 0.11;

    y = options.grid.yp;
    u     = 1 + (Delta_U/2)*tanh(2*y) + ...
            eps1*(1-tanh(y/2).^2).*cos(n1*y)*sin(om1*t) + ...
            eps2*(1-tanh(y/2).^2).*cos(n2*y)*sin(om2*t);


else
    u = zeros(length(x)*length(y),1);
end

end