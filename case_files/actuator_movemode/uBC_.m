function u = uBC_(x,y,t,options)
% boundary conditions for u for actuator


%%

% coordinate left side domain:
x1 = options.grid.x1;

y1 = options.grid.y1;
y2 = options.grid.y2;

y_m = (y1+y2)/2;
y_l = y2-y1;

t_end = options.time.t_end;

% Y1 = @(t) y1  + y_l/2 * t/t_end;
% Y2 = @(t) y_m + y_l/2 * t/t_end;

Y1 = @(t) y1  - y_l + 2*y_l * t/t_end;
Y2 = @(t) y1        + 2*y_l * t/t_end;

u  = zeros(length(x)*length(y),1);

% inflow:
if (length(x)==1 && abs(x-x1)<1e-10)

    u = .1*(y>=Y1(t)).*(y<=Y2(t)).*((y-Y1(t)).*(Y2(t)-y));

end


end