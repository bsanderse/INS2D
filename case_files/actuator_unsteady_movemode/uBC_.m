function u = uBC_(x,y,t,options,dare_vectors)
% boundary conditions for u for actuator

% dare_vectors: 1 if you dare to use vectors as inputs, otherwise
if (nargin<5)
    dare_vectors = 0;
end

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

move_parab = @(y) .1*(y>=Y1(t)).*(y<=Y2(t)).*((y-Y1(t)).*(Y2(t)-y));

% inflow:
if (length(x)==1 && abs(x-x1)<1e-10) 

%     u = .1*(y>=Y1(t)).*(y<=Y2(t)).*((y-Y1(t)).*(Y2(t)-y));
u = move_parab(y);

elseif dare_vectors
    if length(x) ~= length(y)
        error("x and y have different length")
    else
        u  = zeros(length(x),1);

        BCys = y(abs(x-x1)<1e-10);
        u(abs(x-x1)<1e-10) = move_parab(BCys);
    end

end


end