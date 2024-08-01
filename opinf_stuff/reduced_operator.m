function X_r = reduced_operator(X)

% assume X has dimension N x N hat

N = size(X,1);

Diff = X(:,1:N);
Conv = X(:,N+1:end);

Conv_r = reduced_convection_operator(Conv);

X_r = [Diff Conv_r];