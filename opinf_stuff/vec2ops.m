function [D,C] = vec2ops(vec,M)

% num_el = @(r) r*(r+r^2);

% options = optimoptions('Display', 'off');
% M = round(fmincon(@(r) norm(num_el(r) - numel(vec)),1),options);

O = reshape(vec,M+M^2,M)';

D = O(:,1:M);
C = O(:,M+1:end);