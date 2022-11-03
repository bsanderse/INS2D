function [conv_diff, Jac_conv_diff_T, Jac_conv_diff_V] = conv_diff_temperature(T,V,t,options,getJacobian)
% evaluate convection and diffusion terms for temperature equation and optionally Jacobian
% the Jacobian is with respect to both T and V
% note that all terms are in integrated (finite volume) form

% visc = options.case.visc;

[conv, Jac_conv_T, Jac_conv_V] = convection_temperature(T,V,t,options,getJacobian);
[diff, Jac_diff_T] = diffusion_temperature(T,t,options,getJacobian);

conv_diff = -conv + diff;

        
Jac_conv_diff_T = -Jac_conv_T + Jac_diff_T;

Jac_conv_diff_V = -Jac_conv_V;
    

end

