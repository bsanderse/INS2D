function [M1_inds,anti_M1_inds,M1_inds_1,anti_M1_inds_1] = incremental_inds(M1,M)
% given the opinf operator for dimension M, computes:
% M1_inds - the indices of the entries of the opinf operator 
%               for dimension M1
% anti_M1-inds - all other indices
% M1_inds_1 - column indices is M1 entries in M operator
% anti_M1_inds_1 - all other column indices

OM1 = zeros(M,M,1+M);

is = 1:M1;
OM1(is,is,1:M1+1) = 1;

OM1 = reshape(OM1,M,M+M^2);
anti_OM1 = ones(M,M+M^2) - OM1;

inds = vec2op(1:M^2*(1+M),M);

M1_inds = OM1.*inds;
anti_M1_inds = anti_OM1.*inds;

M1_inds = nonzeros(M1_inds');
anti_M1_inds = nonzeros(anti_M1_inds');

M1_inds_1 = M1_inds(M1_inds <= M+M^2);
anti_M1_inds_1 = anti_M1_inds(anti_M1_inds <= M+M^2);





