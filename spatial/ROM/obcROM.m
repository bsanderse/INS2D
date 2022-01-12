function [y_O_ROM, Jac_y_O_ROM] = obcROM(R,t,options,getJacobian)

Jac_y_O_ROM = -666;
if getJacobian
    error('Jacobian not implemented in offline decomposition')
end

R_inhom = get_a_inhom(t,options);

y_O_ROM =   options.rom.obc_hom*kron(R,R)   ...
          + options.rom.obc_inhom*kron(R_inhom,R_inhom) ...
%           + options.rom.obc_hom_inhom*kron(R,R_inhom)   ...
%           + options.rom.obc_inhom_hom*kron(R_inhom,R);
          + options.rom.obc_hom_inhom2*kron(R,R_inhom);
            
% test -> if nonzero, fix convection!
%  norm(options.rom.obc_hom_inhom*kron(R,R_inhom)   ...
%      + options.rom.obc_inhom_hom*kron(R_inhom,R) ...
%      - options.rom.obc_hom_inhom2*kron(R,R_inhom))
% explanation: in obc_hom_inhom2, the submatrices of obc_inhom_hom are not
% included as connected submatrices but kronecker product-like spread over
% the submatrices of obc_hom_inhom. This way, we differ from the
% description in the Master thesis