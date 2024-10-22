function [Conv, Jac] = convectionROM_split(R_convecting,R_convected,t,options,getJacobian)
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

        if norm(Conv_linear) ~= 0
            error("implementation of reduced convection operator not ready to deal with differing convecting and convected velocities")
        end
        
        Conv        = Conv_quad*kron(R_convecting,R_convected) + Conv_linear*R_convecting + yConv;
         
        if (getJacobian == 1)          
           % d/dR (kron(R,R)) = kron(E,R) + kron(R,E), where E = speye(M)
           E   = speye(M);
           % Jac = Conv_quad*(kron(E,R_convected) + kron(R_convecting,E)) + Conv_linear;
           Jac = kron(R_convecting,E)) + Conv_linear;
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