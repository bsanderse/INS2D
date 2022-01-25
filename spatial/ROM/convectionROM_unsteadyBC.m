function [Conv, Jac] = convectionROM_unsteadyBC(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

if getJacobian
    error('Jacobian not implemented in offline decomposition')
else
    Jac = -666;
end

R_inhom = options.rom.abc(t);
R_bc    = [options.rom.abc1(t); options.rom.abc2(t)];

Conv = options.rom.C_hom*kron(R,R) ...
     + options.rom.C_hom_inhom*kron(R,R_inhom) ...
     + options.rom.C_hom_bc*kron(R,R_bc) ...
     + options.rom.C_inhom*kron(R_inhom,R_inhom) ...
     + options.rom.C_inhom_bc*kron(R_inhom,R_bc) ...
     + options.rom.C_bc*kron(R_bc,R_bc);