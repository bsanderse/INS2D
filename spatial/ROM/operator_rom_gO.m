function [gO_hom,gO_inhom] = operator_rom_gO(P,options)

factor = options.BC.gO(0);

B = options.rom.B;
phi_inhom = options.rom.phi_inhom;

id_normal = options.grid.id_normal;
id_tangential = options.grid.id_tangential;
id_n_t = id_normal+id_tangential;

gO_hom = factor*P*(id_n_t.*B);
gO_inhom = factor*P*(id_n_t.*phi_inhom);