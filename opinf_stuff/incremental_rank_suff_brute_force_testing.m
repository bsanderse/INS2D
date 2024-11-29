clear all

% r = options.rom.M;
r = 8;

basis0 = eye(r);
basis = basis0;
for s = 1:r
    offset = zeros(r,1);
    offset(s) = 1; 
    basis = [basis, basis0 + offset]; % note: this basis includes redundant kronecker ambiguity pairs
end

columns_s = [];
for s = 1:r
    columns = [s s*r+(1:s)];
    columns_s = [columns_s columns];
end

basis2 = basis(:,columns_s);
