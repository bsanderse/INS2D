function f = EC_stimulating_function(t,V)

r = size(V,1);

r_hat = r*(r-1)/2;
primes_ = primes(r_hat^2);
primes_ = primes_(1:r_hat);

B_sum = zeros(r,r);
for s = 1:r
    for t = s+1:r
        B = spalloc(r,r,2);
        B(s,t) = 1;
        B(t,s) = -1;
        
        B_sum = B_sum + sin(t*primes_(r*(s-1)+t))*B;
    end
end

f = B_sum*V;