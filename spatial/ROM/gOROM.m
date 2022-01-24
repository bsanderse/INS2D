function gO_ROM = gOROM(R,t,options)

R_inhom = get_a_inhom(t,options);

gO_ROM = options.rom.gO_hom*R  ...
       + options.rom.gO_inhom*R_inhom;
   
end