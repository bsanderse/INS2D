function O = vec2op(vec,M)

O = reshape(vec,M+M^2,M)';
