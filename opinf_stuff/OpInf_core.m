% perform core of OpInf: solve least squares problem and decompose solution
% into linear and quadratic operator
function [L,Q] = OpInf_core(Uhats,RHSs)

% Uhat = [Us' vectorwise_kron(Us)'];

O = (Uhats'\RHSs')';

r = size(RHSs,1);

L = O(:,1:r);
Q = O(:,r+1:end);