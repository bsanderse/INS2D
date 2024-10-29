function f = stimulating_function(t,r)

B0 = eye(r);
B = B0;
for s = 1:r
    offset = zeros(r,1);
    % offset(s) = 1;  % apparently, this offset is not needed to stimulate full rank
    B = [B, B0 + offset];
end

% B
r_hat = size(B,2);
primes_ = primes(r_hat^2);
primes_ = primes_(1:r_hat);

f = B*sin(t*primes_)';
% f = B*(t.^(1:r_hat))';   % not good
% f = B*cos(t./primes_)';  % also not good
