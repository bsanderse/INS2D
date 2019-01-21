%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy
figure
plot(time,k/k(1),'s-')
grid

%% vorticity
% compare e.g. with PhD thesis figure 3.4

figure
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx+1,Npy+1);
% for Re=1000: labels = -4:0.5:4;
labels= 20;
contour(x,y,omega',labels);
axis square
colorbar
grid