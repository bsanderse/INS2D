function [gO_ROM, Jac_gO_ROM] = gOROM(R,t,options,getJacobian)

if getJacobian    
    Jac_gO_ROM = options.rom.gO_hom;
else
    Jac_gO_ROM = -666;
end

gO_ROM = options.rom.gO_hom*R;

if options.rom.rom_bc == 2 || options.rom.rom_bc == 1
    if options.rom.bc_recon == 3
        R_inhom = get_a_inhom(t,options);
        
        gO_ROM = gO_ROM  ...
            + options.rom.gO_inhom*R_inhom;
    elseif options.rom.bc_recon == 2 || options.rom.bc_recon == 0 || ...
            options.rom.bc_recon == 4 || options.rom.bc_recon == 5
        % nothing to do
    else
        error('Sorry, obc offline decomposition not implemented for chosen rom_bc')
    end
% else
%             % nothing to do
end



   
end