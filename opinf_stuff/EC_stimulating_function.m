function f = EC_stimulating_function(t,V)

r = size(V,1);

r_hat = r*(r-1)/2;
primes_ = primes(10*r^2);
primes_ = primes_(1:r_hat);

B_sum = zeros(r,r);
counter = 1;
for s = 1:r
    for t = s+1:r
        B = spalloc(r,r,2);
        B(s,t) = 1;
        B(t,s) = -1;
        
        B_sum = B_sum + sin(t*primes_(counter))*B;
        counter = counter +1;
    end
end

% norm(B_sum-B_sum')

f = B_sum*V;