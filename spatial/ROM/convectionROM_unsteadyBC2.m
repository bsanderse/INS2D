function [Conv, Jac] = convectionROM_unsteadyBC2(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians
Jac = -666;
if getJacobian
    error('Jacobian not implemented in offline decomposition')
end

R_inhom = get_a_inhom(t,options);
R_bc    = get_a_bc(t,options);

Conv = options.rom.C_hom2*kron(R,R) ...
     + options.rom.C_hom_inhom2*kron(R,R_inhom) ...
     + options.rom.C_hom_bc2*kron(R,R_bc) ...
     + options.rom.C_inhom2*kron(R_inhom,R_inhom) ...
     + options.rom.C_inhom_bc2*kron(R_inhom,R_bc) ...
     + options.rom.C_bc2*kron(R_bc,R_bc);