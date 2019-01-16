% load from Xfoil data
addpath('force/airfoil_data');
airfoil_Cp = load('cp_naca0012_a0.dat');
airfoil_Cf = load('cf_naca0012_a0.dat');
airfoil_coord = load('naca0012.dat');

np   = 80; % number of panels on one side

Cp   = squeeze(airfoil_Cp(:,2));
xCp  = squeeze(airfoil_Cp(:,1));
x_k  = squeeze(airfoil_coord(:,1));
y_k  = squeeze(airfoil_coord(:,2));
Cf   = squeeze(airfoil_Cf(:,2));
xCf  = squeeze(airfoil_Cf(:,1));

if (max(abs(xCp(1:np)-xCp(2*np:-1:np+1)))>1e-12)
    error('x-distribution upper and lower not equal');
else
    Cpu  = Cp(1:np);
    Cpl  = Cp(np+1:2*np);
    xCp  = xCp(1:np);
end

if (max(abs(xCf(1:np)-xCf(np+1:2*np)))>1e-12)
    error('x-distribution upper and lower not equal');
else
    Cfu  = Cf(1:np);
    Cfl  = Cf(np+1:2*np);
    xCf  = xCf(1:np);
end


% determine for each airfoil point the surrounding points and weights
xn = zeros(np,2);
yn = zeros(np,2);
for i=1:np
    [xn(i,:) yn(i,:)] = nearest_neighbours(x_k(i),y_k(i),xin,yp);

end
