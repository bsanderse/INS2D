function [d2, Jac] = diffusionROM_unsteadyBC(R,t,options,getJacobian)

visc = options.case.visc;

M = options.rom.M;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;

abc = options.rom.abc;

Diff  = options.rom.Diff;
DiffBC  = options.rom.DiffBC;
yDiff = options.rom.yDiff;

Jac = spalloc(M,M,0);

vbc1 = get_unsteadyVbc(t,options);
vbc2 = options.rom.Vbc0*options.rom.abc(t);
norm(vbc1-vbc2)

switch visc
    
    case 'laminar'
        d2     = Diff*R + DiffBC*abc(t) + yDiff;
%         d2     = Diff*R; %pfusch
        
        if (getJacobian == 1)
            Jac    = Diff;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

end

