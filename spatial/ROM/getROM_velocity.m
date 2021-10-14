function [R] = getROM_velocity(V,t,options)
% get ROM velocity coefficients

% subtract boundary condition contribution (zero if not used)
% if V is a NV*Nt matrix, then this vector is subtracted from each column
V   = V - options.rom.Vbc;

switch options.rom.rom_type
    
    case {'POD', 'Fourier'} % orthogonal basis
        
        B   = options.rom.B;

        if (options.rom.weighted_norm == 0)
            R   = B'*V;
        elseif (options.rom.weighted_norm == 1)
            Om = options.grid.Om;
            R  = B'*(Om.*V);
        end
        
    case 'FDG' % non-orthogonal basis, oblique projection
        NV = options.grid.NV;

%         if (options.rom.weighted_norm == 0)
%             Diag = ones(NV,1);
%         elseif (options.rom.weighted_norm == 1)
%             Diag = options.grid.Om_inv;
%         end
%         Om = spdiags(Diag,NV,NV);
        R  = options.rom.B_inv \ (options.rom.B' * (options.grid.Om .* V));
%         % laplacian
%         L = options.rom.C'*Om*options.rom.C;
%         % we have nabla^2 psi = -omega,
%         % where omega = curl V
%         R = options.rom.Phi' * (L \ (options.rom.C'*(Om*V)));
        
    otherwise
        
        error('wrong choice for ROM type');
        
end

