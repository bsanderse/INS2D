function a_inhom = get_a_inhom(t,options)

    R_inhom = options.rom.R_inhom;
    
    a_inhom = R_inhom*get_a_bc(t,options);