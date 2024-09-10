function vec = symm_neg_def_diffusion(vec,M)
% takes the whole (diffusion+convection) vectorized operator and returns a
% whole vectorized operator with symmetric positive definite diffusion

[D,C] = vec2ops(vec,M);
D2 = symm_neg_def_diffusion_matrix(D);
vec = ops2vec(D2,C);

