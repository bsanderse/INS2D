%% post-processing shear layer
xu = options.grid.xu;
yu = options.grid.yu;
xv = options.grid.xv;
yv = options.grid.yv;
xpp = options.grid.xpp;
ypp = options.grid.ypp;

%% kinetic energy
figure
plot(time,k/k(1),'s-')
grid

%% vorticity
figure
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx,Npy);
contour(xp,yp,omega');
