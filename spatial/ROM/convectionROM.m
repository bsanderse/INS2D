function [Conv, Jac] = convectionROM(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

order4     = options.discretization.order4;
regularize = options.case.regularize;

M   = options.rom.M;
Jac = spalloc(M,M,0);



if (order4==0)
    
    %% no regularization
    if (regularize == 0)
        
        Conv_quad   = options.rom.Conv_quad;
        Conv_linear = options.rom.Conv_linear;
        yConv       = options.rom.yConv;
        
        Conv        = Conv_quad*kron(R,R) + Conv_linear*R + yConv;
         
    else
        error('not implemented');
    end
 

else
    
    error('not implemented');
    
end