function [val,ind] = min2d(A)

% finds minimum and indices for a 2-dimensional matrix A
% can be easily extended to more dimensions

[val, minindex] = min(A(:));    % Find the minimum element in A
                                % The minimum value is MinValue, the index is MinIndex
[ind1, ind2] = ind2sub(size(A), minindex); % Convert MinIndex to subscripts
ind = [ind1,ind2];