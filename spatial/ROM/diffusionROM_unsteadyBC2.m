function [d2, Jac] = diffusionROM_unsteadyBC2(R,t,options,getJacobian)

visc = options.case.visc;

Diff  = options.rom.Diff2;

switch visc
    
    case 'laminar'
        d2     = Diff*R;
        
        if options.rom.rom_bc == 2 && options.rom.bc_recon == 3
            R_inhom = get_a_inhom(t,options);
            R_bc    = get_a_bc(t,options);
            
            DiffBC  = options.rom.DiffBC2;
            yDiff = options.rom.yDiff2;
            
            d2 = d2 + DiffBC*R_inhom + yDiff*R_bc;           
        end
        
        if (getJacobian == 1)
            Jac    = Diff;
        else
            Jac = -666;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

end

