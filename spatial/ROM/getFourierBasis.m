function [Phi] = getFourierBasis(options,M)
% GETBASIS: get ROM spatial basis
% input: mesh; number of modes
% output: spatial basis Phi, 

% set-up discrete Fourier transforms

% assume a signal (vector) u_i, i=1..N is given (possibly complex), 
% the DFT used here is given by
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
M  = options.rom.M;

%%
[~,phi_real_x_trunc_inv] = DFT_matrix_real(Nx,sqrt(M));
[~,phi_real_y_trunc_inv] = DFT_matrix_real(Ny,sqrt(M));

% Phi as used in the code is actually the inverse truncated real DFT, i.e. 
% u \approx Phi*u_hat
Phi     = kron(phi_real_y_trunc_inv,phi_real_x_trunc_inv);
% note that Phi_inv, as given by
% Phi_inv = kron(phi_real_y_trunc,phi_real_x_trunc);
% is constructed such that Phi'*Phi_inv = I_M
