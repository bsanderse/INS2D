function [Conv, Jac] = convectionROM(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

order4     = options.discretization.order4;
regularize = options.case.regularize;

M   = options.rom.M;



if (order4==0)
    
    %% no regularization
    if (regularize == 0)
        
        Conv_quad   = options.rom.Conv_quad;
        Conv_linear = options.rom.Conv_linear;
        yConv       = options.rom.yConv;
        
        Conv        = Conv_quad*kron(R,R) + Conv_linear*R + yConv;
         
        if (getJacobian == 1)          
           % d/dR (kron(R,R)) = kron(E,R) + kron(R,E), where E = speye(M)
           E   = speye(M);
           Jac = Conv_quad*(kron(E,R) + kron(R,E)) + Conv_linear;
        else
           Jac = 0;
%            Jac = spalloc(M,M,0);
        end
        
    else
        error('not implemented');
    end
 

else
    
    error('not implemented');
    
end