function [z,dz] = nonuniform_grid(deltaz,z_low,z_up,sz)
% generate a non-uniform grid, from z_low to z_up, starting deltaz and
% having stretch factor close to sz
% v2.0: includes check for uniform grid and, in that case, adapts deltaz if
% necessary

eps  = 1e-12;

if (sz==1)
    n = (z_up-z_low)/deltaz;
    if (abs(n-round(n))>eps)
        deltaz = (z_up-z_low)/floor(n);
    end
    z  = (z_low:deltaz:z_up)';   
    dz = diff(z);
    return;
end
    
i    = 1;
z(i) = z_low;

% first check the number of grid points by using the specified sz
while (z(i)+deltaz*(sz^i))<=z_up
    i    = i+1;
    z(i) = z(i-1) + deltaz*(sz^(i-2)); 
end

% adapt sz such that the upper boundary is exactly satisfied
% sum of powerseries should be z_up-z_low
S = z_up-z(1);           % sum wanted
n = length(z)-1;         % number of intervals

%Secant method for nonlinear iteration
s(1) = sz;
s(2) = sz^2;
i    = 1;
while (abs(s(i+1)-s(i))>eps)
%     s(i+2) = s(i+1) - (s(i+1)-s(i))*(S*(s(i+1)-1) - deltaz*s(i+1)*(s(i+1)^n-1)) / ...
%                       ( (S*(s(i+1)-1) - deltaz*s(i+1)*(s(i+1)^n-1)) - (S*(s(i)-1) - deltaz*s(i)*(s(i)^n-1)));
    s(i+2) = s(i+1) - (s(i+1)-s(i))*( S*(s(i+1)-1) - deltaz*(s(i+1)^n-1) ) / ...
                      ( (S*(s(i+1)-1) - deltaz*(s(i+1)^n-1)) - (S*(s(i)-1) - deltaz*(s(i)^n-1)));


    i=i+1;
end
sz = s(end);

%check
% deltaz*(sz^n-1)/(sz-1);s


%update z now based on the new sz
z      = zeros(n,1);
i      = 1;
z(i,1) = z_low;
while (i<n+1)
    z(i+1,1)  = z(i)+deltaz*(sz^(i-1)); 
    i         = i+1;
end

dz   = diff(z);

z  = z(:);
dz = dz(:);