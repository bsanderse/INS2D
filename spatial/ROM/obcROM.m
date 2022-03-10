function [y_O_ROM, Jac_y_O_ROM] = obcROM(R,t,options,getJacobian)

if getJacobian
    M   = options.rom.M;
    E   = speye(M);
    
    Jac_y_O_ROM = options.rom.obc_hom*kron(E,R) + options.rom.obc_hom*kron(R,E);
    
    if options.rom.rom_bc == 2 || options.rom.rom_bc == 1
        if options.rom.bc_recon == 3
            R_inhom = get_a_inhom(t,options);
            
            Jac_y_O_ROM = Jac_y_O_ROM + options.rom.obc_hom_inhom2*kron(E,R_inhom);
        elseif options.rom.bc_recon == 2
            % nothing to do
        else
            error('Sorry, precomputation not implemented for bc_recon =/= 3')
        end
    end
else
    Jac_y_O_ROM = -666;
end

y_O_ROM =   options.rom.obc_hom*kron(R,R);

if options.rom.rom_bc == 2 || options.rom.rom_bc == 1
    if options.rom.bc_recon == 3
        
        R_inhom = get_a_inhom(t,options);
        
        y_O_ROM =   y_O_ROM   ...
            + options.rom.obc_inhom*kron(R_inhom,R_inhom) ...
            ... %           + options.rom.obc_hom_inhom*kron(R,R_inhom)   ...
            ... %           + options.rom.obc_inhom_hom*kron(R_inhom,R);
            + options.rom.obc_hom_inhom2*kron(R,R_inhom);
    elseif options.rom.bc_recon == 2
        % nothing to do
    else
        error('Sorry, precomputation not implemented for bc_recon =/= 3')
    end
end

%%
% test -> if nonzero, fix convection!
%  norm(options.rom.obc_hom_inhom*kron(R,R_inhom)   ...
%      + options.rom.obc_inhom_hom*kron(R_inhom,R) ...
%      - options.rom.obc_hom_inhom2*kron(R,R_inhom))
% explanation: in obc_hom_inhom2, the submatrices of obc_inhom_hom are not
% included as connected submatrices but kronecker product-like spread over
% the submatrices of obc_hom_inhom. This way, we differ from the
% description in the Master thesis