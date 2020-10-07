function [d2, Jac] = diffusionROM(R,t,options,getJacobian)

visc = options.case.visc;

M = options.rom.M;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;

Diff  = options.rom.Diff;
yDiff = options.rom.yDiff;

Jac = spalloc(M,M,0);

switch visc
    
    case 'laminar'
        d2     = Diff*R + yDiff;
        
        if (getJacobian == 1)
            Jac    = Diff;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

end

