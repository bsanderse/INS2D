function [x,y] = actuator_unsteady_mesh(options)

deltax = options.grid.deltax;
sx     = options.grid.sx;
x1     = options.grid.x1;
x2     = options.grid.x2;

deltay = options.grid.deltay;
sy     = options.grid.sy;
y1     = options.grid.y1;
y2     = options.grid.y2;

x      = nonuniform_grid(deltax,x1,x2,sx);
y      = nonuniform_grid(deltay,y1,y2,sy);

% transform uniform grid to non-uniform cosine grid
% L_x    = x2-x1;
% L_y    = y2-y1;
% x      = (L_x/2) * (1 - cos(pi*x/L_x));
% y      = (L_y/2) * (1 - cos(pi*y/L_y));

end