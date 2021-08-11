function a_bc = get_a_bc(t,options)

    phi_bc = options.rom.phi_bc;
    
%     a_bc = phi_bc'*get_bc_vector_yBC(t,options);
      a_bc = phi_bc'*get_bc_vector_yBC(t,options); %pfusch