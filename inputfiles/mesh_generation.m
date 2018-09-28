%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generation of a (non-)uniform mesh                      

% x and y give the position of the grid lines, i.e. where the velocities are defined
% 
    x      = nonuniform_grid(deltax,x1,x2,sx);
%     x      = nonuniform_grid2(x1,x2,sx,Nx);

%     x_1    = nonuniform_grid2(x1,x2/2,(2/pi)*Nx,Nx/2);
%     x_2    = nonuniform_grid2(x2/2,x2,10,Nx/2);
%     x      = [x_1(1:end-1); flipud(1-x_1)];

%     x  = nonuniform_grid(2/Re,x1,x2,1.05);
%         x_1 = nonuniform_grid(deltax,x1,x_c-0.2,1.05);
%         x_1 = x_c-0.2 + flipud(-x_1);
%         x_2 = nonuniform_grid(deltax,x_c-0.2,x_c+1.2,1);
%         x_3 = nonuniform_grid(deltax,x_c+1.2,x2,1.05);
%         x = [x_1; x_2(2:end); x_3(2:end)];
    %     x_1           = flipud(x_c-x_1);
%         [x_1,dx]      = nonuniform_grid(deltax,0,x2,1.01);
%         x = [flipud(-x_1); x_1(2:end)];
    %     x             = [x_1(1:end-1); x_2; x_3(2:end)];

    y      = nonuniform_grid(deltay,y1,y2,sy);
%     y      = nonuniform_grid2(y1,y2,sy,Ny);

%     y = x;
%         [y_1,~]      = nonuniform_grid(deltay,y_c,y_c+0.15,1);
%         [y_2,~]      = nonuniform_grid(deltay,y_c+0.15,y2,1.05);
%         y_3 = [y_1; y_2(2:end)];
%         y   = [y_c-flipud(y_3-y_c); y_3(2:end)];

%         [y_2,dy]      = nonuniform_grid(deltay,1,L_y/2,1.05);

    %     y             = [y_1; y_2(2:end)];
    %     y             = [-flipud(y); y(2:end)];
    %     y             = y+L_y/2;
    %     [y_2,dy]      = nonuniform_grid(deltay,1/2,L_y/2,sy);
    %     y_3           = [y_1; y_2(2:end)];
    %     y             = [-flipud(y_3(2:end));y_3];

    %     y             = [-flipud(y_1(2:end)); y_1];

    %     y = [y1 y2(2:end)];
    %     y = [fliplr(-y_1) y_1(2:end)];

    % y = y+L_y/2;
    % transform uniform grid to non-uniform cosine grid
    x = (L_x/2) * (1 - cos(pi*x/L_x));
    y = (L_y/2) * (1 - cos(pi*y/L_y));
%     x = L_x * (1 - cos(0.5*pi*x/L_x));
%     y = L_y * (1 - cos(0.5*pi*y/L_y));
    % y = y-L_y/2;
%     x = sin(pi*x/L_x);
%     y = sin(pi*y/L_y);
    
%     x_2 = nonuniform_grid(deltax,x2,30,1.005);
%     x = [x;x_2(2:end)];

    % refine each cell in both directions
    % x = refine(x);
    % y = refine(y);
    % coarsen each cell in both directions
    % x = coarsen(x);
    % y = coarsen(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% mesh quantities

% pressure positions
xp          = (x(1:end-1)+x(2:end))/2;
yp          = (y(1:end-1)+y(2:end))/2;

% distance between velocity points
hx          = diff(x);
hy          = diff(y);

Nx          = length(hx); 
Ny          = length(hy);

% distance between pressure points
gx          = zeros(Nx+1,1);
gx(1)       = hx(1)/2;
gx(2:Nx)    = (hx(1:Nx-1)+hx(2:Nx))/2;
gx(Nx+1)    = hx(end)/2;
 
gy          = zeros(Ny+1,1);
gy(1)       = hy(1)/2;
gy(2:Ny)    = (hy(1:Ny-1)+hy(2:Ny))/2;
gy(Ny+1)    = hy(end)/2;


fprintf(fcw,['Nx=' num2str(Nx) ', min(hx)=' num2str(min(hx)) '\n']); 
fprintf(fcw,['Ny=' num2str(Ny) ', min(hy)=' num2str(min(hy)) '\n']);
