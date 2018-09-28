% prescribing Cp along chord

% load from Xfoil data
% airfoil_Cp = load('cp_naca0012_a0.dat');
% airfoil_Cf = load('cf_naca0012_a0.dat');
% airfoil_coord = load('naca0012.dat');
% 
% Cp   = squeeze(airfoil_Cp(:,2));
% xCp  = squeeze(airfoil_Cp(:,1));
% x_k  = squeeze(airfoil_coord(:,1));
% y_k  = squeeze(airfoil_coord(:,2));
% Cf   = squeeze(airfoil_Cf(:,2));
% xCf  = squeeze(airfoil_Cf(:,1));
% 
% if (max(abs(xCp(1:80)-xCp(160:-1:81)))>1e-12)
%     error('x-distribution upper and lower not equal');
% else
%     dCp  = Cp(1:80)-Cp(160:-1:81);
%     xCp  = xCp(1:80);
% end
% 
% if (max(abs(xCf(1:80)-xCf(81:160)))>1e-12)
%     error('x-distribution upper and lower not equal');
% else
%     Cfu  = Cf(1:80);
%     Cfl  = Cf(81:160);
%     xCf  = xCf(1:80);
% end

% own panel method
nk        = 1000; 
nk2       = floor(nk/2);
i         = 1:nk2+1;
x_kb      = sin(0.5*pi*(i-1)*(1/nk2)).^2;   
x_ko      = sin(0.5*pi*(1.0-(i-1.0)*(1/nk2))).^2;  

% naca 0012
y_kb      = +0.6*(0.29690*x_kb.^0.5 - 0.12600*x_kb - 0.35160*x_kb.^2 + 0.28430*x_kb.^3 -0.10360*x_kb.^4);
y_ko      = -0.6*(0.29690*x_ko.^0.5 - 0.12600*x_ko - 0.35160*x_ko.^2 + 0.28430*x_ko.^3 -0.10360*x_ko.^4);

% cylinder
% y_kb      =  sqrt(0.25-(x_kb-0.5).^2);
% y_ko      = -sqrt(0.25-(x_ko-0.5).^2);


x_k       = [x_kb(1:end-1) x_ko(1:end)];
y_k       = [y_kb(1:end-1) y_ko(1:end)];
aoa       = 2;
[Cp,Vs,xCp,yCp] = Panel(x_k,y_k,aoa);  % call panel method to calculate Cp distribution
dCp       = Cp(1:nk2)-Cp(end:-1:nk2+1);  % difference between upper and lower
xCp       = xCp(1:nk2);

npoints    = 1/deltax;
xmesh      = (0:deltax:1)';
cp_mesh    = interp1(xCp,dCp,xmesh,'spline','extrap');

% change cp_mesh so that area under curve remains the same?

ymid       = floor(Nvy_in/2)+1;
y1Dx       = zeros(Nvx_in,1);
x_le       = floor(Nvx_in/3); % leading edge position
y1Dx(x_le:x_le+length(xmesh)-1) = hx(x_le:x_le+length(xmesh)-1).*cp_mesh;
y1Dy       = zeros(Nvy_in,1);
y1Dy(ymid) = 1;

Fy         = kron(y1Dy,y1Dx);

% cfu_mesh   = interp1(xCf,Cfu,xmesh,'spline','extrap');
% cfl_mesh   = interp1(xCf,Cfl,xmesh,'spline','extrap');
% y1Dx_l     = zeros(Nux_in,1);
% y1Dx_u     = zeros(Nux_in,1);
% x_le       = floor(Nux_in/3);
% 
% y1Dx_l(x_le:x_le+length(xmesh)-1) = hx(x_le:x_le+length(xmesh)-1).*cfl_mesh;
% y1Dx_u(x_le:x_le+length(xmesh)-1) = hx(x_le:x_le+length(xmesh)-1).*cfu_mesh;
% 
% y1Dy_l      = zeros(Nuy_in,1);
% y1Dy_u      = zeros(Nuy_in,1);
% y1Dy_l(ymid-1)   = -1;
% y1Dy_u(ymid) = -1;
% Fx         = kron(y1Dy_l,y1Dx_l) + kron(y1Dy_u,y1Dx_u);