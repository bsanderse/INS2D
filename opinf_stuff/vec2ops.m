function [D,C] = vec2ops(vec)

num_el = @(r) r*(r+r^2);

M = round(fmincon(@(r) norm(num_el(r) - numel(vec)),1));

O = reshape(vec,M+M^2,M)';

D = O(:,1:M);
C = O(:,M+1:end);