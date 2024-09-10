function [A_dot,A] = time_difference_quotient(A_raw, method,dt, skip)
% skip: output only every skip-th vector pair of A_dot and A


switch method
    case "forward euler"
        A_dot = (A_raw(:,2:end,:) - A_raw(:,1:end-1,:))/dt;
        A = A_raw(:,1:end-1,:);

end

A_dot = A_dot(:,1:skip:end);
A = A(:,1:skip:end);

