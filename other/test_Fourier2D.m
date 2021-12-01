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

clearvars
close all

%%
Nx = 12; % number of points, assume EVEN
Lx = 2*pi;
Ny = 8; % number of points, assume EVEN
Ly = 2*pi;
x  = linspace(0,Lx-Lx/Nx,Nx)';
y  = linspace(0,Ly-Ly/Ny,Ny)';
% dx = x(2)-x(1);
% dy = y(2)-y(1);
[X,Y] = ndgrid(x,y);

I = sqrt(-1);

% freq = 1/dx; %sampling frequency should be at least 2*highest frequency in signal
% u    = 2 + 0.1*cos(2*pi*5*(X-0.2).^2 + 10*pi/180) + 0.6*sin(2*pi*120*X) + sin(2*pi*Y.^2) ;
u     = sin(X).*sin(Y);
Nx    = length(x);
Ny    = length(y);


%% truncation
M = 2;


%% method 1: use Matlab FFT
% for physical interpretation, one could add scaling of 1/N in the definition of uhat
% if scaling is included it needs to be added also in ifft
% we use the 1/sqrt(N) scaling to get an orthonormal basis
u_hat = fft2(u)/sqrt(Nx*Ny); 

% inverse transform
u_back = ifft2(u_hat)*sqrt(Nx*Ny);

%% method 2: program the DFT matrix Phi and Phi_inv, such that
% u_hat = phi_x*u*phi_y

% and similarly u = phi_x_inv*u_hat*
% phi     = zeros(Nx,Nx); 
% phi_inv = zeros(Nx,Nx);


[phi_x,phi_x_inv] = DFT_matrix(Nx);
[phi_y,phi_y_inv] = DFT_matrix(Ny);

u_hat2  = phi_x*u*phi_y;
u_back2 = phi_x_inv*u_hat2*phi_y_inv;


% alternative
% phi_x * u * phi_y = (phi_x *u) * phi_y = (phi_y.' * (phi_x*u).').'
% = (phi_y * (phi_x*u).').' = DFT( DFT(u).').';
u_hat5  = DFT( DFT(u).').';
u_back5 = IDFT( IDFT(u_hat5.') .');

% full 2D matrix form
% note that this uses the symmetry of phi_y
% phi_full = kron(phi_y,speye(Nx))*kron(speye(Ny),phi_x);
phi_2D   = kron(phi_y,phi_x);
u_hat4   = reshape(phi_2D*u(:),Nx,Ny);
phi_inv_2D =  phi_2D'; % *1/(Nx*Ny)
u_back4  = reshape(phi_inv_2D*u_hat4(:),Nx,Ny);


% note that phi' = N*phi_inv, where ' is the complex conjugate transpose
% furthermore, phi*phi_inv = phi_inv*phi = I
% and phi*phi' = phi'*phi = N*I

% uhat2 should equal uhat:
max(max(abs(u_hat2-u_hat)))
max(max(abs(u_hat5-u_hat)))
max(max(abs(u_hat4-u_hat)))

% uback2 should equal uback
max(max(abs(u_back2-u_back)))
max(max(abs(u_back5-u_back)))
max(max(abs(u_back4-u_back)))

% note the ordering of uhat:
% uhat[1] (k=0): mean of solution
% uhat[2...N/2+1]: positive frequencies
% uhat[N/2+2..N]: negative frequencies
% can use fftshift to change this

% for real data u, we get
% uhat[1]: mean of solution
% uhat[2...N/2+1] = conj(uhat[N/2+2..N])

%% method 3: since we have real data, we can use the real DFT
% the coefficients u_hat will appear as complex conjugates
% note Phi = exp(-I*2*pi*k*i/N) = cos(2*pi*k*i/N) - I*sin(2*pi*k*i/N)
% so 
% u_hat = Phi*u = Phi_Re*u - I*Phi_Im*u
% now u_hat appears in complex conjugates as follows
% u_hat(1): mean, real-valued
% u_hat(2:N/2) : complex 
% u_hat(N/2+1): real
% u_hat(N/2+2:N) : conj of uhat(2:N/2) when put in reverse order, i.e.
% conj(u_hat(end:-1:N/2+2))- u_hat(2:N/2) = 0
% we can therefore just store N/2+1 coefficients so we take only the first
% N/2+1 rows of Phi


% subsequently, we can store the real and imaginary parts of u_hat separately
% we store N/2+2 real part values and N/2-1 imaginary part values

% u_hat = Phi_real*u
% and similarly u = Phi_real_inv*u_hat

[phi_real_x, phi_real_x_inv] = DFT_matrix_real(Nx);
[phi_real_y, phi_real_y_inv] = DFT_matrix_real(Ny);

% phi_full = kron(phi_real_y,speye(Nx))*kron(speye(Ny),phi_real_x);
phi_real_2D = kron(phi_real_y,phi_real_x);

% the resulting phi_real*u gives coefficients u_real that are not the same
% as u_hat, but related as follows: [real(u_hat(1:N/2+1));imag(u_hat(N/2+2:N))];

% alternatively, we can use a matrix that selects the right modes from the
% original matrix Phi to get Phi_real
% phi_real2 = [real(phi_x(1:Nx/2+1,:)); imag(phi_x(Nx/2+2:Nx,:))];
% max(max(abs(phi_real-phi_real2)))

%%
% the inverse transform is then as follows
% u = Phi_inv * uhat

phi_real_inv_2D = kron(phi_real_y_inv,phi_real_x_inv);

u_hat3  = phi_real_2D*u(:);
u_back3 = phi_real_inv_2D*u_hat3;

max(abs(u_back3-u_back(:)))

% alternatively:
u_hat6  = RDFT( RDFT(u).').';
u_back6  = IRDFT( IRDFT(u_hat6.').');


max(abs(u_hat6(:)-u_hat3(:)))
max(abs(u_back6(:)-u_back(:)))


%% truncation of phi and phi_real

% truncate the complex exponential form
% index truncation set, simply based on frequency ordering (not on PSD)

[phi_x_trunc,phi_x_trunc_inv] = DFT_matrix(Nx,M);

% u_hat_trunc   = phi_x_trunc*u;
% u_back_trunc  = phi_x_trunc_inv*u_hat_trunc;

[phi_y_trunc,phi_y_trunc_inv] = DFT_matrix(Ny,M);
phi_trunc_2D  = kron(phi_y_trunc,phi_x_trunc);
phi_trunc_inv_2D  = kron(phi_y_trunc_inv,phi_x_trunc_inv);

u_hat_trunc   = phi_trunc_2D*u(:);
u_back_trunc  = phi_trunc_inv_2D*u_hat_trunc;

u_hat_trunc2  = DFT(DFT(u,M).',M).';
u_hat_trunc2  = DFT2(u,M,M);
u_back_trunc2 = IDFT(IDFT(u_hat_trunc2.',Ny).',Nx);
u_back_trunc2 = IDFT2(u_hat_trunc2,Nx,Ny);

max(abs(u_hat_trunc(:) - u_hat_trunc2(:)))
max(abs(u_back_trunc2(:) - u_back_trunc(:)))


%% truncate the real form
% phi_real_inv_trunc = [phi_real_x_inv(:,1:M+1) phi_real_x_inv(:,end-M+1:end)];
% phi_real_trunc = [phi_real_x(1:M+1,:); phi_real_x(end-M+1:end,:)];
% u_hat_trunc3  = phi_real_trunc*u;
% u_back_trunc3  = phi_real_inv_trunc*u_hat_trunc3;

[phi_real_x_trunc,phi_real_x_trunc_inv] = DFT_matrix_real(Nx,2*M);

% u_hat_trunc2   = phi_real_x_trunc*u;
% u_back_trunc2  = phi_real_x_trunc_inv*u_hat_trunc2;

% we should still have phi*phi_inv = I_M
% max(max(abs(phi_real_trunc*phi_real_inv_trunc - eye(2*M+1))))
% but phi_inv*phi does NOT equal I_N
% phi_real_inv_trunc*phi_real_trunc

% max(abs(u_back_trunc2(:) - u_back_trunc(:)))

%
[phi_real_y_trunc,phi_real_y_trunc_inv] = DFT_matrix_real(Ny,2*M);

phi_real_trunc_2D     = kron(phi_real_y_trunc,phi_real_x_trunc);
phi_real_trunc_inv_2D = kron(phi_real_y_trunc_inv,phi_real_x_trunc_inv);

u_hat7  = phi_real_trunc_2D*u(:);
u_back7 = phi_real_trunc_inv_2D*u_hat7;

u_hat8  = RDFT(RDFT(u,2*M).',2*M).';
u_hat8  = RDFT2(u,2*M,2*M);
u_back8 = IRDFT(IRDFT(u_hat8.',Ny).',Nx);
u_back8 = IRDFT2(u_hat8,Nx,Ny);

max(abs(u_back7 - u_back_trunc(:)))

max(abs(u_hat8(:) - u_hat7(:)))
max(abs(u_back8(:) - u_back7(:)))

%% 
figure
plotstyle = {'LineWidth',1,'MarkerSize',10};
surf(x,y,reshape(u,Nx,Ny)')

hold on
surf(x,y,reshape(u_back5,Nx,Ny)')
% plot(u_back3,'o',plotstyle{:})
% plot(real(u_back_trunc),plotstyle{:})
% plot(u_back_trunc3,'s',plotstyle{:})
xlabel('x')
ylabel('y')
legend('original','real DFT+iDFT truncated');

figure
surf(x,y,u'-reshape(u_back5,Nx,Ny)');
title('error')
xlabel('x')
ylabel('y')

%% plot solution in frequency domain
figure
surf(1:Nx,1:Ny,reshape(u_hat3,Nx,Ny)')
figure
Mx = size(phi_real_x_trunc,1);
My = size(phi_real_y_trunc,1);
surf(1:Mx,1:My,reshape(u_hat7,Mx,My)')


%%
% psd  = u_hat.*conj(u_hat); %/(n^2); % n^2 disappears if it we do fft/n
% 
% % power in time and in frequency domain (should match)
% norm(u)^2/Nx;
% sum(psd);
% 
% fVals = freq*(0:Nx/2-1)/Nx;
% figure
% semilogy(fVals,psd(1:Nx/2))
% 
% % select indices with largest PSD by sorting the PSD
% [val,ind] = sort(psd,'desc');
% %
% ind_select = ones(Nx,1);
% % select 2*n_keep indices, where the factor 2 is needed because we need the
% % coefficients associated with negative frequencies as well to do the inverse
% % transform
% n_keep = 3;
% ind_select(ind(2*n_keep:end)) = 0;
% % set other indices to 0
% fhat_new = ind_select.*u_hat;
% fnew    = Nx*ifft(fhat_new,Nx);
% 
% % note: if we want physical meaning out of the Fourier coefficients
% % (amplitude, angle) we need to divide by n (and multiply by 2 if we don't
% % include the complex conjugates)
% ind_keep = ind(1:2*n_keep-1);
% f_keep = u_hat(ind_keep);
% abs(f_keep);
% 
% % the phase angle is more tricky business
% % https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/
% %angle(f_keep)
% % f2 = f_keep;
% % threshold = max(abs(f_keep))/10000; %tolerance threshold
% % f2(abs(f_keep)<threshold) = 0; %maskout values that are below the threshold
% % phase=atan2(imag(f2),real(f2))*180/pi; %phase information
% 
% % get only positive frequencies and mean
% ind_pos = 2:2:2*(n_keep-1);
% abs(u_hat(1))
% abs(u_hat(ind_keep(ind_pos)))*2
% f2 = u_hat;
% threshold = max(abs(f2))/1e4; %tolerance threshold
% f2(abs(f2)<threshold) = 0; %maskout values that are below the threshold
% angle(f2(ind_keep(ind_pos)))*180/pi
% 
% 
% figure
% plot(x,u);
% hold on
% plot(x,fnew);