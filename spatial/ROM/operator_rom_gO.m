function [gO_hom,gO_inhom] = operator_rom_gO(P,options)

factor = options.BC.gO(0);

B = options.rom.B;
phi_inhom = options.rom.phi_inhom;

gO_hom = factor*P*B;
gO_inhom = factor*P*phi_inhom;