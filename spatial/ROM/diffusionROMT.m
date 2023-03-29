function [d2, Jac] = diffusionROMT(RT,t,options,getJacobian)

visc = options.case.visc;

MT = options.rom.MT;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;

DiffT  = options.rom.DiffT;
yDiffT = options.rom.yDiffT;

Jac = spalloc(MT,MT,0);

switch visc
    
    case 'laminar'
        d2     = DiffT*RT + yDiffT;
        
        if (getJacobian == 1)
            Jac    = Diff;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

% diffphi=reshape(DiffT,NT*NT,1);
% save("ROM_diffphi.dat", "d2", "-ascii");
% disp("Hi")

end

