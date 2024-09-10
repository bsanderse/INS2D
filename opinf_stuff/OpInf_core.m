% perform core of OpInf: solve least squares problem and decompose solution
% into linear and quadratic operator
function [L,Q] = OpInf_core(Uhats,RHSs)

% Uhat = [Us' vectorwise_kron(Us)'];

O = (Uhats'\RHSs')';

r = size(RHSs,1);

L = O(:,1:r);
Q = O(:,r+1:end);

% just to test

% [A,b,T] = matrix_lsq_to_vector_ls(Uhats,RHSs);
% 
% c = A\b;
% 
% O2 = zeros(size(T(:,:,1)));
% for i = 1:numel(O2)
%     O2 = O2 + c(i)*T(:,:,i);
% end
% 
% norm(O-O2)
% 
% r