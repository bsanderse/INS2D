function [RT] = getROM_Temperature(T,t,options)
% get ROM temperature coefficients

%switch options.rom.rom_type
    
 %   case {'POD', 'Fourier'} % orthogonal basis
        
        BT = options.rom.BT;

        if (options.rom.weighted_norm_T == 0)
            RT   = BT'*T;
        elseif (options.rom.weighted_norm_T == 1)
            Om = options.grid.Om;
            RT  = BT'*(Om.*T);
        end
  %  otherwise
        
   %     error('wrong choice for ROM type');
        
%end
    