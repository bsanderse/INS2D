function u_hat_real = RDFT2(u,Nx_hat,Ny_hat,vec_trunc)
% RDFT2: return 2D real discrete Fourier transform
% the definition here is such that the transform is orthonormal
% the input u should be of size Nx*Ny

% the truncation parameter N_hat is used to select a subset of rows from B (cols
% of Binv), obtained by ordering in terms of frequencies.
% this is done in such a way that the resulting matrix B has roughly size N_hat*N
% note that each frequency has a positive and a negative component, which
% are both selected
% therefore N_hat is basically 2*the number of frequencies ('modes') that are selected
% and one needs N_hat to be divisible by 2
% example:
% * N_hat=1 means only k=0 frequency, output size of B = 1*N
% * N_hat=2 means k=0, k=1 and k=N-1 (positive and negative frequencies),
%   output size of B = 3*N
% * etc.

% the truncation vector vec_trunc can be used instead of N_hat to indicate
% which rows are to be selected, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies

if (nargin==1)
    u_hat_real = RDFT(RDFT(u).').';
elseif (nargin==2)
    u_hat_real = RDFT(RDFT(u,Nx_hat).',Nx_hat).';
elseif (nargin==3)
    u_hat_real = RDFT(RDFT(u,Nx_hat).',Ny_hat).';
else
    error('truncation based on vector indices is not finished for 2D');
end
    

end