function [Conv, Jac] = convectionROM_unsteadyBC2(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

if getJacobian    
    M   = options.rom.M;
    E   = speye(M);

    Jac = options.rom.C_hom2*(kron(E,R)+kron(R,E));
else
    Jac = -666;
end

Conv = options.rom.C_hom2*kron(R,R);
    
if options.rom.rom_bc == 2 && options.rom.bc_recon == 3
    R_inhom = get_a_inhom(t,options);
    R_bc    = get_a_bc(t,options);
    
    Conv = Conv ...
        + options.rom.C_hom_inhom2*kron(R,R_inhom) ...
        + options.rom.C_hom_bc2*kron(R,R_bc) ...
        + options.rom.C_inhom2*kron(R_inhom,R_inhom) ...
        + options.rom.C_inhom_bc2*kron(R_inhom,R_bc) ...
        + options.rom.C_bc2*kron(R_bc,R_bc);
    
    if getJacobian
        Jac = Jac ...
            + options.rom.C_hom_inhom2*kron(E,R_inhom) ...
            + options.rom.C_hom_bc2*kron(E,R_bc);
    end
end