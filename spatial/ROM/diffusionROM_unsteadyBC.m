function [d2, Jac] = diffusionROM_unsteadyBC(R,t,options,getJacobian)

visc = options.case.visc;

M = options.rom.M;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;

abc = options.rom.abc;
abc1 = options.rom.abc1;
abc2 = options.rom.abc2;

Diff  = options.rom.Diff;
DiffBC  = options.rom.DiffBC;
yDiff = options.rom.yDiff;

Jac = 0;%spalloc(M,M,0);

% vbc1 = get_unsteadyVbc(t,options);
% vbc2 = options.rom.Vbc0*options.rom.abc(t);
% norm(vbc1-vbc2)

switch visc
    
    case 'laminar'
%         d2     = Diff*R + DiffBC*abc(t) + yDiff*abc(t);
        d2     = Diff*R + DiffBC*abc(t) + yDiff(:,1)*abc1(t) + yDiff(:,2)*abc2(t);
%         d2     = Diff*R; %pfusch
        
        if (getJacobian == 1)
            Jac    = Diff;
        end
        
    otherwise
        error('turbulent diffusion not implemented for ROM');
        
end

% Diffu  = options.discretization.Diffu;
% Diffv  = options.discretization.Diffv;
% yDiffu = options.discretization.yDiffu;
% yDiffv = options.discretization.yDiffv;
% 
% D_h = blkdiag(Diffu, Diffv);
% y_D = [yDiffu; yDiffv];
% 
% B   = options.rom.B;
% 
% V_r = getFOM_velocity(R,t,options);
% norm(B'*(D_h*V_r+y_D)-d2)
% norm(B'*D_h*B-Diff)
% norm(B'*D_h*get_unsteadyVbc(t,options)-DiffBC*abc(t))
% % norm(B'*D_h*y_D-yDiff*abc(t))
% % norm(B'*y_D-(yDiff*abc(t)))
% norm(B'*y_D-(yDiff(:,1)*abc1(t) + yDiff(:,2)*abc2(t)))


end

