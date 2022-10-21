function [Phi] = getFourierBasis(options,M)
% GETFOURIERBASIS: get ROM spatial basis based on INVERSE Fourier transform
% input: mesh; number of modes (optional)
% output: spatial basis Phi, such that u \approx Phi*uhat,
% with uhat the Fourier coefficients

% assume a signal (vector) u_i, i=1..N is given (possibly complex), 
% the DFT is given by
% u_k = sum_i=0^(N-1) u_i exp(-I*2*pi*k*i/N), k=0..N-1
% where u_k are complex valued coefficients
% this definition is consistent with the Matlab definition
% https://www.mathworks.com/help/matlab/ref/fft.html

% alternatively, one can use k=-N/2..N/2-1 (see e.g. Wikipedia)
% https://en.wikipedia.org/wiki/Discrete_Fourier_transform

% the inverse DFT is given by
% u_i = (1/N) sum_k=0^(N-1) u_k exp(-I*2*pi*k*i/N), i=0..N-1


%% basis for u-velocity
% x = options.grid.xin;
Nx = options.grid.Nx;
Ny = options.grid.Ny;
% M  = options.rom.M;

%%
% limit to M modes
if (nargin>1 && M>0)
    Mbasis = 2*round(sqrt(M)/2); % round to nearest even number
else
    Mbasis = 0;
end

[~,phi_real_x_trunc_inv] = DFT_matrix_real(Nx,Mbasis);
[~,phi_real_y_trunc_inv] = DFT_matrix_real(Ny,Mbasis);

% Phi as used in the code is actually the inverse truncated real DFT, i.e. 
% u \approx Phi*u_hat
Phi     = kron(phi_real_y_trunc_inv,phi_real_x_trunc_inv);
% note that Phi_inv, as given by
% Phi_inv = kron(phi_real_y_trunc,phi_real_x_trunc);
% is constructed such that Phi'*Phi_inv = I_M
