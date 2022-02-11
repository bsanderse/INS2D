function [V] = getFOM_velocity(R,t,options)
% get FOM velocity based on ROM coefficients

B   = options.rom.B;
if (options.rom.rom_bc == 2)
    Vbc = get_unsteadyVbc(t,options);
else
    Vbc = options.rom.Vbc(:,1);
end

V   = B*R + Vbc;

