function [d2, Jac] = diffusionROM_unsteadyBC(R,t,options,getJacobian)

visc = options.case.visc;

abc = options.rom.abc;
abc1 = options.rom.abc1;
abc2 = options.rom.abc2;

Diff  = options.rom.Diff;
DiffBC  = options.rom.DiffBC;
yDiff = options.rom.yDiff;



switch visc
    
    case 'laminar'
        d2     = Diff*R + DiffBC*abc(t) + yDiff(:,1)*abc1(t) + yDiff(:,2)*abc2(t);
        
        if (getJacobian == 1)
            Jac    = Diff;
        else
            Jac = -666;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

end

