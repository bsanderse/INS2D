%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy
% figure
% plot(time,(k-k(1))/k(1),'s-')
% grid
cd ../..
open('results/shear_layer_ROM/energy_error_ROM_inviscid_momcons.fig');
k0 = snapshots.k(1);
hold on
semilogy(time,abs(k-k0)/k0,'s--')
saveas(gcf,'results/shear_layer_ROM/energy_error_ROM_inviscid_momcons.fig');


%% momentum
% figure
% plot(time,umom,'s-')
% hold on
% plot(time,vmom,'o-')

open('results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');
hold on
semilogy(time,abs(umom-snapshots.umom(1))/snapshots.umom(1),'s--')
saveas(gcf,'results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');

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