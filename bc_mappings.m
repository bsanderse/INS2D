function options = bc_mappings(options)

yBC = get_bc_vector_yBC(0,options);

Nbc = numel(yBC);

options.grid.Nbc = Nbc;
yBC0 = zeros(Nbc,1);
yM0 = get_yM(options,yBC0);

%% F_M mapping yBC to yM
F_M = zeros(options.grid.Np,Nbc);

for i =1:Nbc
    yBC1 = yBC0;
    yBC1(i) = 1;
    yM1 = get_yM(options,yBC1);
    F_M(:,i) = yM1 - yM0; % (divide by Delta x = 1)
end
options.discretization.F_M = sparse(F_M);

%% similarly, one could compute the other mappings