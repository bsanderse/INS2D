function [x,y] = shear_layer_mesh(options)

Nx     = options.grid.Nx;
sx     = options.grid.sx;
x1     = options.grid.x1;
x2     = options.grid.x2;
L_x    = x2-x1;
deltax = L_x/Nx;               % uniform mesh size x-direction                                   

Ny     = options.grid.Ny;
sy     = options.grid.sy;
y1     = options.grid.y1;
y2     = options.grid.y2;
L_y    = y2-y1;
deltay = L_y/Ny;               % uniform mesh size y-direction

x      = nonuniform_grid(deltax,x1,x2,sx);
y      = nonuniform_grid(deltay,y1,y2,sy);

% transform uniform grid to non-uniform cosine grid
% x      = (L_x/2) * (1 - cos(pi*x/L_x));
% y      = (L_y/2) * (1 - cos(pi*y/L_y));

end