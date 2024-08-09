function O_ = incremental_opinf(A_hat,A_dot_T,operator_constraint, constraint_rhs, increments)
% increments: list of ints (in increasing order) smaller than M 
%               for which sub-operators get computed

% current design choice: gets called AFTER constraints are defined
% to avoid recomputing constraints

M = size(A_dot_T,2);

%% August 2024 notation
A = A_hat';
B = A_dot_T;
X = zeros(M+M^2,M);
%%

O_ = zeros(M*(M+M^2),1);

% not for loop, but recursive calls!
% for M1 = increments
% recursion root
if isempty(increments)
    [operator_constraint,constraint_rhs] = get_constraints(M,"EC-OpInf skew"); % botch!
    O_ = lsqlin(sparse(kron(eye(M),(A_hat'))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
    O_
else
    M1 = increments(end);
    [M1_inds,anti_M1_inds,M1_inds_1,anti_M1_inds_1] = incremental_inds(M1,M);

    A_hat_M1 = A_hat(M1_inds_1,:);
    A_dot_T_M1 = A_dot_T(:,1:M1);
    operator_constraint_M1 = operator_constraint(:,M1_inds);
    % constraint_rhs_M1 = constraint_rhs(M1_inds,:);
    constraint_rhs_M1 = constraint_rhs;

    O_M1 = incremental_opinf(A_hat_M1,A_dot_T_M1, ...
        operator_constraint_M1, constraint_rhs_M1, increments(1:end-1));

    O_(M1_inds) = O_M1;

    X = vec2op(O_,M)';
    B_M1 = A*X; % contribution of M1-sub-operator

    A_inc1 = sparse(kron(eye(M1),A(:,anti_M1_inds_1)));
    A_inc2 = sparse(kron(M-M1,A));
    A_inc = blkdiag(A_inc1,A_inc2);
    B = B - B_M1;
    B_vec = B(:);
    % B_inc = B_vec(anti_M1_inds_1);
    B_inc = B_vec;
    operator_constraint_inc = operator_constraint(:,anti_M1_inds);
    % constraint_rhs_inc = constraint_rhs(:,anti_M1_inds);
    constraint_rhs_inc = constraint_rhs;

    O_inc_vec = lsqlin(A_inc,B_inc, [],[], operator_constraint_inc, constraint_rhs_inc);

    O_(anti_M1_inds) = O_inc_vec;
end



