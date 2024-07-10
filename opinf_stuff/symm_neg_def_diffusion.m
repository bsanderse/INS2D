function vec = symm_neg_def_diffusion(vec,M)
% takes the whole (diffusion+convection) vectorized operator and returns a
% whole vectorized operator with symmetric psoitive definite diffusion

vec_D = vec(1:M^2);
D = reshape(vec_D,M,M);
new_D = -D*D';

vec(1:M^2) = new_D(:);



%% test
% D = diag([1 2 3])
% D2 = symm_neg_def_diffusion(D(:),3)
% norm(D2 - diag(-[1 4 9]));