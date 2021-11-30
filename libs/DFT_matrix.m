function [B,Binv] = DFT_matrix(N,Ntrunc,vec_trunc)
% DFT_matrix: return 1D discrete Fourier transform basis and its inverse
% the definition here is such that
% B*B_inv = B_inv*B = I

% the truncation parameter Ntrunc is used to select a subset of rows from B (cols
% of Binv), obtained by ordering in terms of frequencies.
% this is done in such a way that the resulting matrix B has roughly size Ntrunc*N
% note that each frequency has a positive and a negative component, which
% are both selected
% therefore Ntrunc is basically 2*the number of frequencies ('modes') that are selected
% and one needs Ntrunc to be divisible by 2
% example:
% * Ntrunc=1 means only k=0 frequency, output size of B = 1*N
% * Ntrunc=2 means k=0, k=1 and k=N-1 (positive and negative frequencies),
%   output size of B = 3*N
% * etc.

% the truncation vector vec_trunc can be used instead of Ntrunc to indicate
% which rows are to be selected, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies


k = 0:N-1; % frequency array
j = 0:N-1; % spatial index

[K,J] = ndgrid(k,j);
I     = sqrt(-1);

B    = (1/sqrt(N))*exp(-I*2*pi*(K.*J)/N);
Binv = (1/sqrt(N))*exp(I*2*pi*(K.*J)/N);
% alternative (equivalent): Binv = B'; with ' Hermitian transpose

% truncate the transform
if (nargin>1)
    
    if (nargin==2)
    % truncate to Ntrunc modes

        if (Ntrunc<=0)
            error(['check value for truncation parameter, Ntrunc=' Ntrunc]);
        elseif (Ntrunc==N)
            warning('truncation parameter equal to N, no truncation');
            ind_trunc = 1:N;
        else
            % check whether Ntrunc is odd, but leave the case of Ntrunc=1 intact
            if (Ntrunc~=1 && mod(Ntrunc,2) ~= 0)
                Ntrunc = 2*ceil(Ntrunc/2);
                warning(['truncation parameter is not even, changing to Ntrunc=' Ntrunc]);
            end
            
            % single mode: ind = 1
            % two modes: ind = [1:2, N]
            % three modes: ind = [1:3, N-1:N]
            ind_trunc = [1:Ntrunc/2+1 (N-Ntrunc/2+1):N];
        end
            
        
    elseif (nargin==3)
    % truncate according to vec_trunc

        ind_trunc = vec_trunc;
    
        
    end
    
    B         = B(ind_trunc,:);
    Binv      = Binv(:,ind_trunc);
            
    
end