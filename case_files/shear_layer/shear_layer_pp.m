%% post-processing shear layer
xu = options.grid.xu;
yu = options.grid.yu;
xv = options.grid.xv;
yv = options.grid.yv;
xpp = options.grid.xpp;
ypp = options.grid.ypp;


figure
plot(time,k/k(1),'s-')
grid