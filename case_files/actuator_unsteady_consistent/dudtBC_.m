function dudt = dudtBC_(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;


if ( (length(x)==1 && abs(x-x1)<1e-10)) %|| (length(y)==1 && abs(y-y1)<1e-10) )
    f     = 0.5;
    a     = pi/6;
    %     alpha = (pi/6)*sin(f*t); % original
    %     u     = cos(alpha)*ones(length(x)*length(y),1); %original
    alpha = a*sin(y-f*t);
%     u     = cos(alpha).*ones(length(x)*length(y),1);
    dudt  = a*(-f)*cos(y-f*t).*(-sin(alpha)).*ones(length(x)*length(y),1);
else
    dudt  = zeros(length(x)*length(y),1);
end


end



