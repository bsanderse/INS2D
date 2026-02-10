function phi = feature(basis)

[r,~] = size(basis);

r_star = r_star_lin_quad(r);

phi = zeros(r_star);

ext_basis = [ones(1,r_star); basis]; % extend basis by constant row;

row = 1;
for r_ = 1:r
    for i = 1:r_+1
        phi(row,:) = basis(r_,:).*ext_basis(i,:);
        row = row +1;
    end
end



