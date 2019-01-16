%% generation of a (non-)uniform mesh                      

file_name = [options.case.project '_mesh'];

if (exist(file_name,'file'))
    
    func  = str2func(file_name);    
    [x,y] = func(options);

else
    
    error(['mesh file ' file_name ' not available']);
    
end


%% derived mesh quantities

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


options.grid.x  = x;
options.grid.y  = y;
options.grid.xp = xp;
options.grid.yp = yp;
options.grid.hx = hx;
options.grid.hy = hy;
options.grid.gx = gx;
options.grid.gy = gy;


fprintf(fcw,['Nx=' num2str(Nx) ', min(hx)=' num2str(min(hx)) '\n']); 
fprintf(fcw,['Ny=' num2str(Ny) ', min(hy)=' num2str(min(hy)) '\n']);