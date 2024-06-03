% function [energies_gauss,energies_linim] = inviscid_energy_sims(Conv)

nt = 10;
dt = 100*options.time.dt;

energies_gauss = zeros(nt,n_rom_types,2);
energies_linim = zeros(nt,n_rom_types,2);
energies_euler = zeros(nt,n_rom_types,2);

% figure(17)

for i = 1:n_rom_types
    for j =1:2
        Conv = Convs(:,:,i,j);

        [energies_gauss(:,i,j),X_gauss] = inviscid_energy_sim(Conv,dt,a0s(:,1),nt,"Gauss method");
        [energies_linim(:,i,j),X_linim] = inviscid_energy_sim(Conv,dt,a0s(:,1),nt,"linear implicit");
        [energies_euler(:,i,j),X_euler] = inviscid_energy_sim(Conv,dt,a0s(:,1),nt,"forward euler");

        % figure(71)
        % plot(energies_linim(:,i,j),'x-','DisplayName',labels_combined(i,j) + " linear implicit")
        % hold on
        figure(72)
        plot(energies_gauss(:,i,j),'+-','DisplayName',labels_combined(i,j) + " Gauss method")
        hold on
        % plot(energies_euler(:,i,j),'+-','DisplayName',labels_combined(i,j) + " forward euler")
        % hold on
    end
end

legend("show")
title("kinetic energy conservation")