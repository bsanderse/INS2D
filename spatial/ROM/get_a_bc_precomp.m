function a_bc = get_a_bc_precomp(t,options)

    phi_bc = options.rom.phi_bc;
    
    a_bc = phi_bc'*get_bc_vector_yBC(t,options);