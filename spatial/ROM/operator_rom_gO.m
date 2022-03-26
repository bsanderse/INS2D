function [gO_hom,gO_inhom] = operator_rom_gO(P,options)

factor = options.BC.gO(0);

B = options.rom.B;

% id_normal = options.grid.id_normal;
% id_tangential = options.grid.id_tangential;
% id_n_t = id_normal+id_tangential;

gO_factor = options.grid.gO_factor;
NV = options.grid.NV;


gO_hom = factor*P*spdiags(gO_factor,0,NV,NV)*B;


if options.rom.rom_bc == 2 || options.rom.rom_bc == 1
    if options.rom.bc_recon == 3
        phi_inhom = options.rom.phi_inhom;
        gO_inhom = factor*P*spdiags(gO_factor,0,NV,NV)*phi_inhom;
    elseif options.rom.bc_recon == 2 || options.rom.bc_recon == 0 || ...
            options.rom.bc_recon == 4 || options.rom.bc_recon == 5
        gO_inhom = -666;
    else
        error('Sorry, obc offline decomposition not implemented for rom_bc == 2, but bc_recon =/=3')
    end
else
    gO_inhom = -666;
end