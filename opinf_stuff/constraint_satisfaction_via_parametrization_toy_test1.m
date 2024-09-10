% toy test: constraint satisfaction via parametrization
% 
% We want to fit a 3 x 3 matrix M to data with the constraint that M is
% skew-symmetric

% basis of skew-symmetric 3 x 3 matrices
P = 3;

T1 = [0 1 0; -1 0 0; 0 0 0];
T2 = [0 0 1; 0 0 0; -1 0 0];
T3 = [0 0 0; 0 0 1; 0 -1 0];
T(:,:,1) = T1;
T(:,:,2) = T2;
T(:,:,3) = T3;

% true M 
ctrue = [3 7 9]';
Mtrue = ctrue(1)*T1 + ctrue(2)*T2 + ctrue(3)*T3;

% data
K = 3;

D = eye(3);
d1 = D(:,1);
d2 = D(:,2);
d3 = D(:,3);

% x_dots
X_dot = Mtrue*D;
xd1 = Mtrue*d1;
xd2 = Mtrue*d2;
xd3 = Mtrue*d3;

% construct linear system
A = zeros(P);
b = zeros(P,1);
for p=1:P
    for k=1:K
        b(p) = b(p) + X_dot(:,k)'*T(:,:,p)*D(:,k);
    end
    for s=1:P
        for k=1:K
            A(p,s) = A(p,s) + D(:,k)'*T(:,:,s)'*T(:,:,p)*D(:,k);
        end
    end
end

c = A\b



