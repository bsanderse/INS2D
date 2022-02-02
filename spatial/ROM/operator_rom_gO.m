function [gO_hom,gO_inhom] = operator_rom_gO(P,options)

factor = options.BC.gO(0);

B = options.rom.B;

id_normal = options.grid.id_normal;
id_tangential = options.grid.id_tangential;
id_n_t = id_normal+id_tangential;

gO_hom = factor*P*diag(id_n_t)*B;


if options.rom.rom_bc == 2
    if options.rom.bc_recon == 3
phi_inhom = options.rom.phi_inhom;
gO_inhom = factor*P*diag(id_n_t)*phi_inhom;
    else
        error('Sorry, obc offline decomposition not implemented for rom_bc == 2, but bc_recon =/=3')
    end
else
    gO_inhom = -666;
end