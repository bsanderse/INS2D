function [d2, Jac] = diffusionROM_unsteadyBC2(R,t,options,getJacobian)

visc = options.case.visc;

R_inhom = get_a_inhom(t,options);
R_bc    = get_a_bc(t,options);

Diff  = options.rom.Diff2;
DiffBC  = options.rom.DiffBC2;
yDiff = options.rom.yDiff2;

switch visc
    
    case 'laminar'
        d2     = Diff*R + DiffBC*R_inhom + yDiff*R_bc;
        
        if (getJacobian == 1)
            Jac    = Diff;
        else
            Jac = -666;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

end

