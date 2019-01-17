psi = get_streamfunction(V,t,options);

x1 = options.grid.x1;
x2 = options.grid.x2;
y1 = options.grid.y1;
y2 = options.grid.y2;

figure(4)
% Re=1000:
labels = [1.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5 1e-6 0 -1e-10 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
% Re=5000:
% labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
% Re=10000:
% labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];

% BFS, Re=800
% labels = [-0.03 -0.025 -0.02 -0.015 -0.01 -0.005 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.49 0.5 0.502 0.504];
% labels = 60;
% labels =-0.12:0.005:0;
% labels1 = 1-(0:0.025:1).^2;
% labels2 = 1+(0:0.025:1).^2;
% labels = [flipud(labels1(2:end)');labels2'];
% labels=0.45:0.002:0.54;
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels);

axis equal
axis([x1 x2 y1 y2]);

xlabel('x');
ylabel('y');

% actuator disks
% hold on
% contour(xu,yu,reshape(Fx,Nux_in,Nuy_in),min(Fx),'Color','black');
% h1=plot([x(xmid+1) x(xmid+1)],[yp(yrange(1)) yp(yrange(length(yrange)))],'k-');
% set(h1,'LineWidth',1.5);
% h2=plot([x(2*xmid+1) x(2*xmid+1)],[yp(yrange(1)+floor(nd/2)) yp(yrange(end)+floor(nd/2))],'k-');
% set(h2,'LineWidth',1.5);

% point sources
% hold on
% theta=0:0.01:2*pi;
% plot(x1+L_x/2+deltax*cos(theta),y1+L_y/2+deltax*sin(theta),'k-');
% 
% hold off
