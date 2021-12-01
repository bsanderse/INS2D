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
N  = 20; % number of points, assume EVEN
L = 1;
x  = linspace(0,L-L/N,N)';
dx = x(2)-x(1);

I = sqrt(-1);

freq = 1/dx; %sampling frequency should be at least 2*highest frequency in signal
u    = zeros(N,2);
u(:,1)   = 2 + 0.1*cos(2*pi*5*(x-0.2).^2 + 10*pi/180) + 0.6*sin(2*pi*120*x);
u(:,2)   = 3 + 0.1*cos(2*pi*4*(x-0.7).^2 + 10*pi/180) + 0.6*sin(2*pi*200*x);

N_vec = size(u,2);
N    = length(x);

%% method 1: use Matlab FFT
% for physical interpretation, one could add scaling of 1/N in the definition of uhat
% if scaling is included it needs to be added also in ifft
% note that the Matlab definition uses 1/N in the ifft and 1 in the fft
% to get to the symmetric form with 1/sqrt(N) in both, we need to add
% 1/sqrt(N) in the fft and *sqrt(N) in the ifft
% u_hat = fft(u,N)/sqrt(N); 
u_hat = DFT(u);

% inverse transform
% u_back = ifft(u_hat,N)*sqrt(N);
u_back = IDFT(u_hat);

%% method 2: program the DFT matrix Phi and Phi_inv, such that
% u_hat = Phi*u
% and similarly u = Phi_inv*u_hat

[phi,phi_inv] = DFT_matrix(N);

% for physical interpretation, one could add 1/N in the definition of uhat
% if scaling is included it needs to be added also in ifft

u_hat2  = phi*u;
u_back2 = phi_inv*u_hat2;


% note that phi' = phi_inv, where ' is the complex conjugate transpose
% this is because we use 1/sqrt(N) scaling in definition of the DFT
% furthermore, phi*phi_inv = phi_inv*phi = I

% check orthogonality
max(max(abs(phi*phi_inv-speye(N))))
% check transpose = inverse
max(max(abs(phi'-phi_inv)))


% uhat2 should equal uhat:
max(abs(u_hat2-u_hat))
% uback2 should equal uback
max(abs(u_back2-u_back))

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

[phi_real,phi_real_inv] = DFT_matrix_real(N);
% check orthogonality
max(max(abs(phi_real*phi_real_inv-speye(N))))
% check transpose = inverse
max(max(abs(phi_real.'-phi_real_inv)))

% phi_real = zeros(N,N);
% 
% i  = 0:N-1; % spatial index is from 0..N-1
% k1 = 0:N/2; % frequency array for real part
% k2 = N/2+1:N-1; % frequency array for im part
% 
% % loop over spatial index
% for j = 1:length(i)
%     phi_real(k1+1,j) = cos(2*pi*k1*i(j)/N);
%     phi_real(k2+1,j) = -sin(2*pi*k2*i(j)/N);
% end

% the resulting phi_real*u gives coefficients u_real that are not the same
% as u_hat, but related as follows: [real(u_hat(1:N/2+1));imag(u_hat(N/2+2:N))];

% alternatively, we can use a matrix that selects the right modes from the
% original matrix Phi to get Phi_real
% phi_real2 = [real(phi(1:N/2+1,:)); imag(phi(N/2+2:N,:))];
% max(max(abs(phi_real-phi_real2)))

%%
% the inverse transform is then as follows
% u = Phi_inv * uhat
% phi_real_inv = zeros(N,N);
% 
% i  = 0:N-1; % spatial index array
% k1 = 0:N/2; % frequency array for real part
% k2 = N/2+1:N-1; % frequency array for im part
% 
% % loop over frequencies
% for j = 1:length(k1)
%     phi_real_inv(:,j)     = (2/N)*cos(2*pi*k1(j)*i/N);
% end
% % note! adapt k=0 and k=N/2:
% phi_real_inv(:,1)     = phi_real_inv(:,1)/2;
% phi_real_inv(:,N/2+1) = phi_real_inv(:,N/2+1)/2;
% 
% for j = 1:length(k2)   
%     phi_real_inv(:,j+k2(1)) = -(2/N)*sin(2*pi*k2(j)*i/N);
% end

u_hat3  = phi_real*u;
u_hat33 = RDFT(u);
u_back3 = phi_real_inv*u_hat3;
u_back33 = IRDFT(u_hat33);

max(abs(u_back3-u_back))


%% test differentiation

% see how fourier transforming affects the derivative
% assume we have a signal in Fourier space, which is transformed to
% physical space, then we take derivative, then we transform back to
% Fourier space

% d^2/dx^2
d = ones(N,1);
D = spdiags([d -2*d d],[-1 0 1],N,N);
D(1,end)=1;
D(end,1)=1;

% effective diagonal matrix (diagonalisation)
Dk_complex  = phi*D*phi_inv;
% same result with real transform, but the order can be different
Dk_real     = phi_real*D*phi_real_inv;
Dk_analytic = 2*cos(2*pi*(0:N-1)/N)-2;

% since D is circulant, the diagonalized form can be obtained directly as follows
D_circ      = D(:,1);
eig_complex = sqrt(N)*(phi')*D_circ; % this should correspond to Dk_complex and to eig(D)
% get back the original matrix using eigenvectors and eigenvalue matrix
Dtest_complex = sqrt(N)*(phi')*diag((phi')*D_circ)*phi;
% similarly, for the real case
% however, we don't have a good expression for the eigenvalues in terms of
% phi_real yet, so we use phi for now; however the ordering is not yet correct
% Dtest_real    = sqrt(N)*(phi_real.')*diag((phi')*D_circ)*phi_real;


figure
plot(diag(Dk_real),'x-')
hold on
plot(real(diag(Dk_complex)),'s-');
plot(Dk_analytic,'o-')
legend('real transform','complex transform','analytic');


%% test truncation of phi and phi_real
Ntrunc = 8; % should be even; note that the truncated signal has typically  length Ntrunc+1

% truncate the complex exponential form;
% u_hat_trunc will just be a subset of u_hat, which should match
% u_hat(ind_trunc)
% ind_trunc = [1:Ntrunc/2 (N-Ntrunc/2+1):N];


% index truncation set, simply based on frequency ordering (not on PSD)
% ind_trunc     = [1:M+1 (N-M+1):N]; 
% phi_inv_trunc = phi_inv(:,ind_trunc);
% phi_trunc     = phi(ind_trunc,:);
[phi_trunc,phi_trunc_inv] = DFT_matrix(N,Ntrunc);
u_hat_trunc  = phi_trunc*u;
% or:
u_hat_trunc  = DFT(u,Ntrunc);
% back:
u_back_trunc = phi_trunc_inv*u_hat_trunc;
% or: 
u_back_trunc = IDFT(u_hat_trunc,N);
% the following quantity should be 0:
M = size(phi_trunc,1);
max(max(abs(phi_trunc*phi_trunc_inv - eye(M))))
% the following quantity is generally nonzero, as an error is made due to
% truncation
max(abs(u_back_trunc - u))


% truncate the real form
% phi_real_inv_trunc = [phi_real_inv(:,1:M+1) phi_real_inv(:,end-M+1:end)];
% phi_real_trunc = [phi_real(1:M+1,:); phi_real(end-M+1:end,:)];
[phi_real_trunc,phi_real_trunc_inv] = DFT_matrix_real(N,Ntrunc);
u_hat_trunc3   = phi_real_trunc*u;
u_back_trunc3  = phi_real_trunc_inv*u_hat_trunc3;
% or 
u_hat_trunc3 = RDFT(u,Ntrunc);
% or:
u_back_trunc3 = IRDFT(u_hat_trunc3,N); 

% we should still have phi*phi_inv = I_M
max(max(abs(phi_real_trunc*phi_real_trunc_inv - eye(M))))
% but phi_inv*phi does NOT equal I_N
% phi_real_inv_trunc*phi_real_trunc
max(abs(u_back_trunc3 - u))
% compare with complex exponential truncation:
max(abs(u_back_trunc3 - u_back_trunc))


%% 
figure
title('frequency truncation')
plotstyle ={'LineWidth',1,'MarkerSize',10};
plot(u,plotstyle{:})
hold on
plot(u_back,'s',plotstyle{:})
plot(u_back3,'o',plotstyle{:})
plot(real(u_back_trunc),plotstyle{:})
plot(u_back_trunc3,'s',plotstyle{:})

legend('original','complex DFT+iDFT','real DFT+iDFT','complex DFT+iDFT truncated','real DFT+iDFT truncated');


%% more complex example, where we filter based on PSD, instead of purely on frequencies
% select 2*n_keep indices, where the factor 2 is needed because we need the
% coefficients associated with negative frequencies as well to do the inverse
n_keep = 5;

% note that u_hat is given by fft(u)*sqrt(N)
psd  = u_hat.*conj(u_hat); %/(n^2); % n^2 disappears if it we do fft/n

% power in time and in frequency domain (should match)
norm(u)^2/N;
sum(psd);

fVals = freq*(0:N/2-1)/N;
figure
semilogy(fVals,psd(1:N/2))

% select indices with largest PSD by sorting the PSD
[val,ind]  = sort(psd(:,1),'desc');
ind_select = ones(N,1);
% select 2*n_keep indices, where the factor 2 is needed because we need the
% coefficients associated with negative frequencies as well to do the inverse
% transform
% for now we use the first vector in u to select, and leave the option to
% the index selection for each vector separately for later
ind_select(ind(2*n_keep:end)) = 0;

% set other indices to 0
fhat_new = ind_select.*u_hat;
fnew    = sqrt(N)*ifft(fhat_new,N); % we use sqrt(N) because this has been used in the forward transform as well

% note: if we want physical meaning out of the Fourier coefficients
% (amplitude, angle) we need to divide by n (and multiply by 2 if we don't
% include the complex conjugates)
ind_keep = ind(1:2*n_keep-1);
f_keep = u_hat(ind_keep);
abs(f_keep);

% the phase angle is more tricky business
% https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/
%angle(f_keep)
% f2 = f_keep;
% threshold = max(abs(f_keep))/10000; %tolerance threshold
% f2(abs(f_keep)<threshold) = 0; %maskout values that are below the threshold
% phase=atan2(imag(f2),real(f2))*180/pi; %phase information

% get only positive frequencies and mean
ind_pos = 2:2:2*(n_keep-1);
abs(u_hat(1,:));
abs(u_hat(ind_keep(ind_pos)))*2;
f2 = u_hat;
threshold = max(abs(f2))/1e4; %tolerance threshold
f2(abs(f2)<threshold) = 0; %maskout values that are below the threshold
angle(f2(ind_keep(ind_pos)))*180/pi;

% alternative approach uses the truncated fft matrix, either in complex or
% real form
% complex form:
[phi_trunc_psd,phi_trunc_psd_inv] = DFT_matrix(N,[],ind_keep);
u_hat_trunc_psd   = phi_trunc_psd*u;
u_trunc_psd = phi_trunc_psd_inv*u_hat_trunc_psd;
% alternatives:
% u_hat_trunc_psd = DFT(u,[],ind_keep);
% u_trunc_psd = IDFT(u_hat_trunc_psd,[],ind_keep);

% real form:
% only use the positive frequencies
% we select the cosine and sine terms to be kept after the fft has been
% applied by using ind_keep_real
ind_keep_real = [ind_keep(ind_keep<=N/2); ind_keep(ind_keep<=N/2 & ind_keep>1) + N/2];
[phi_trunc_real_psd,phi_trunc_real_psd_inv] = DFT_matrix_real(N,[],ind_keep_real);
u_hat_real_trunc_psd   = phi_trunc_real_psd*u;
u_real_trunc_psd = phi_trunc_real_psd_inv*u_hat_real_trunc_psd;
% alternatives:
u_hat_real_trunc_psd = RDFT(u,[],ind_keep_real);
u_real_trunc_psd = IRDFT(u_hat_real_trunc_psd,N,ind_keep_real);

figure
title('psd truncation');
plot(x,u);
hold on
plot(x,fnew,'x-');
plot(x,u_trunc_psd,'s-')
plot(x,u_real_trunc_psd,'o-')
legend('original signal','first few modes based on taking fft and then truncating',...
    'first few modes based on using truncated complex fft matrix',...
    'first few modes based on using truncated real fft matrix');