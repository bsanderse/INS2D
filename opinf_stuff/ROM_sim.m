function [A_ROM] = ROM_sim(D,C,a0,dt,Nt)

r = numel(a0);
A_ROM = zeros(r,Nt);
A_ROM(:,1) = a0;
a = a0;
for i = 2:Nt
    kron_a = kron(a,a);
    if size(C,2) < r^2 % handle reduced convection operator
        [~,~,u] = reduced_coordinates(r);
        kron_a = kron_a(u);
    end

    a = a+dt*(D*a+C*kron_a);
    A_ROM(:,i) = a;
end