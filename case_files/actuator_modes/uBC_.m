function u = uBC_(x,y,t,options)
% boundary conditions for u for actuator

% coordinate left side domain:
x1 = options.grid.x1;
% coordinate lower side domain:
% y1 = options.grid.y1;


if ( (length(x)==1 && abs(x-x1)<1e-10)) %|| (length(y)==1 && abs(y-y1)<1e-10) )
    f     = 0.5;
    %     alpha = (pi/6)*sin(f*t); % original
    alpha = (pi/6)*sin(y-f*t);
    %     u     = cos(alpha)*ones(length(x)*length(y),1); %original
    u     = cos(alpha).*ones(length(x)*length(y),1);
else
    u = zeros(length(x)*length(y),1);
end

%%

% coordinate left side domain:
x1 = options.grid.x1;

y1 = options.grid.y1;
y2 = options.grid.y2;


%     x1      = 0;
%     x2      = 15;
%     y1      = -0.5;
%     y2      = 0.5;

y_25 = .25*(y2-y1)+y1;
y_5  = .5 *(y2-y1)+y1;
y_75 = .75*(y2-y1)+y1;

u  = zeros(length(x)*length(y),1);

% inflow:
if (length(x)==1 && abs(x-x1)<1e-10)
%     u = (y>=0).*(24*y.*(1/2-y));

    u1 = (y>=y1).*(y<=y_25).*(24*(y-y1).*(y_25-y));

    u2 = (y>=y_25).*(y<=y_5).*(24*(y-y_25).*(y_5-y));

    u3 = (y>=y_5).*(y<=y_75).*(24*(y-y_5).*(y_75-y));

    u4 = (y>=y_75).*(y<=y2).*(24*(y-y_75).*(y2-y));

    t_end = options.time.t_end;

    time_fac = @(n,t) sin(n*2*pi*t/t_end);
    time_fac = @(n,t) (time_fac(n,t)>=0)*time_fac(n,t)+1;

    u = time_fac(1,t)*u1 + time_fac(2,t)*u2 ...
        + time_fac(3,t)*u3 + time_fac(4,t)*u4;
    u = u/10;
end


end