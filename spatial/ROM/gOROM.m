function [gO_ROM, Jac_gO_ROM] = gOROM(R,t,options,getJacobian)



if getJacobian    
    Jac_gO_ROM = options.rom.gO_hom;
else
    Jac_gO_ROM = -666;
end

if options.rom.rom_bc == 2
    if options.rom.bc_recon == 3
        R_inhom = get_a_inhom(t,options);
    else
        error('Sorry, obc offline decomposition not implemented for rom_bc == 2, but bc_recon =/=3')
    end
else
    R_inhom = 0;
end


gO_ROM = options.rom.gO_hom*R  ...
       + options.rom.gO_inhom*R_inhom;
   
end