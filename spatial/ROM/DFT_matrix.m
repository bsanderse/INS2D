function [B,Binv] = DFT_matrix(N,Ntrunc)
% DFT_matrix: return 1D discrete Fourier transform basis and its inverse
% the definition here is such that
% B*B_inv = B_inv*B = I

k = 0:N-1; % frequency array
j = 0:N-1; % spatial index

[K,J] = ndgrid(k,j);
I     = sqrt(-1);

B    = (1/sqrt(N))*exp(-I*2*pi*(K.*J)/N);
Binv = (1/sqrt(N))*exp(I*2*pi*(K.*J)/N);
% alternative: Binv = B'; with ' Hermitian transpose

% truncate to Ntrunc modes
if (nargin>1 && Ntrunc>0)
    
    ind_trunc = [1:Ntrunc+1 (N-Ntrunc+1):N];
    B         = B(ind_trunc,:);
    Binv      = Binv(:,ind_trunc);    

end