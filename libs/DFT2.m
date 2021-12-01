function u_hat = DFT2(u,Nx_hat,Ny_hat,vec_trunc)
% DFT2: return 2D discrete Fourier transform 
% the definition here is such that the transform is orthonormal

% the input u should be of size Nx*Ny

if (nargin==1)
    u_hat  = DFT(DFT(u).').';
elseif (nargin==2)
    u_hat  = DFT(DFT(u,Nx_hat).',Nx_hat).';
elseif (nargin==3)
    u_hat  = DFT(DFT(u,Nx_hat).',Ny_hat).';
else
    error('truncation based on vector indices is not finished for 2D');
end
    

    
end