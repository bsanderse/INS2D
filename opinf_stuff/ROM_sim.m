function [A_ROM] = ROM_sim(D,C,a0,dt,Nt)

r = numel(a0);
A_ROM = zeros(r,Nt);
A_ROM(:,1) = a0;
a = a0;
for i = 2:Nt
    a = a+dt*(D*a+C*kron(a,a));
    A_ROM(:,i) = a;
end