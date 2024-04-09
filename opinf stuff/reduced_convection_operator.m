function C = reduced_convection_operator(C_unreduced)

N = size(C_unreduced,1);

if size(C_unreduced,2) < N^2 % C_unreduced is actually already reduced
    C = C_unreduced;
else

    s_ = reduced_coordinates(N);

    Nr = max(s_);

    C = zeros(N,Nr);

    for j = 1:N^2
        C(:,s_(j)) = C(:,s_(j)) + C_unreduced(:,j);
    end

end
