function [A_dot,A] = time_difference_quotient(A_raw, method,options)

switch method
    case "forward euler"
        A_dot = (A_raw(:,2:end,:) - A_raw(:,1:end-1,:))/options.time.dt;
        A = A_raw(:,1:end-1,:);

end
