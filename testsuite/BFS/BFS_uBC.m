function u = BFS_uBC(x,y,t,options)
% boundary conditions for u for BFS

% coordinate left side domain:
x1 = options.grid.x1;

u  = zeros(length(x)*length(y),1);

% inflow:
if (length(x)==1 && abs(x-x1)<1e-10)
    u = (y>=0).*(24*y.*(1/2-y));
end

end