%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy
figure
% plot(time,(k-k(1))/k(1),'s-')
% grid
cd ../..
% open('results/shear_layer_ROM/energy_error_ROM_inviscid_momcons.fig');
if (options.rom.rom == 1 )
    k0 = snapshots.k(1);
else
    k0 = k(1);
end
% hold on
semilogy(time,abs(k-k0)/k0,'s--')
title('energy error');
% saveas(gcf,'results/shear_layer_ROM/energy_error_ROM_inviscid_momcons.fig');


%% momentum
figure
% plot(time,umom,'s-')
% hold on
% plot(time,vmom,'o-')
if (options.rom.rom == 1 )
    umom0 = snapshots.umom(1);
else
    umom0 = umom(1);
end
% open('results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');
hold on
semilogy(time,abs(umom-umom0)/umom0,'s--')
title('u-momentum error');
% saveas(gcf,'results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');

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