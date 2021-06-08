function [B,Binv] = DFT_matrix_real(N,Ntrunc)
% DFT_MATRIX_REAL: return 1D discrete Fourier transform basis for real data
% vectors, transforming data vector into a set of (N/2+1) real and (N/2-1) imaginary
% components
% the definition here is such that
% B*B_inv = B_inv*B = I
% 

j  = 0:N-1; % spatial index is from 0..N-1
k1 = 0:N/2; % frequency array for real part
k2 = N/2+1:N-1; % frequency array for im part

[K1,J1] = ndgrid(k1,j);
[K2,J2] = ndgrid(k2,j);

% involve sqrt(2/N) to make Binv = B'
B     = sqrt(2/N)*[cos(2*pi*(K1.*J1)/N); -sin(2*pi*(K2.*J2)/N)];

% note! adapt k=0 and k=N/2:
B(1,:)     = B(1,:)/sqrt(2);
B(N/2+1,:) = B(N/2+1,:)/sqrt(2);

% full expression:
% Binv  = sqrt(2/N)*[cos(2*pi*(K1.*J1).'/N) -sin(2*pi*(K2.*J2).'/N)];
% faster alternative:
Binv  = B'; % note that B is real, so ' and .' are the same

% note! adapt k=0 and k=N/2:
% Binv(:,1)     = Binv(:,1)/sqrt(2);
% Binv(:,N/2+1) = Binv(:,N/2+1)/sqrt(2);

% truncate to Ntrunc/2 modes
if (nargin>1 && Ntrunc>0)
    Ntrunc    = Ntrunc/2;
    B         = [B(1:Ntrunc+1,:); B(end-Ntrunc+1:end,:)];
    Binv      = [Binv(:,1:Ntrunc+1) Binv(:,end-Ntrunc+1:end)];    

end