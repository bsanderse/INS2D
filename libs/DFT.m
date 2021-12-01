function u_hat = DFT(u,N_hat,vec_trunc)
% DFT: return 1D discrete Fourier transform 
% the definition here is such that the transform is orthonormal

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

N     = length(u);
u_hat = fft(u)/sqrt(N);

% truncate the transformed value u_hat
if (nargin>1)
    
    if (nargin==2)
    % truncate to N_hat modes

        if (N_hat<=0)
            error(['check value for truncation parameter, N_hat=' N_hat]);
        elseif (N_hat>=N)
            warning('truncation parameter greater than or equal to N, no truncation');
            ind_trunc = 1:N;
        else
            % check whether N_hat is odd, but leave the case of N_hat=1 intact
            if (N_hat~=1 && mod(N_hat,2) ~= 0)
                N_hat = 2*ceil(N_hat/2);
                warning(['truncation parameter is not even, changing to N_hat=' N_hat]);
            end
            
            % single mode: ind = 1
            % two modes: ind = [1:2, N]
            % three modes: ind = [1:3, N-1:N]
            ind_trunc = [1:N_hat/2+1 (N-N_hat/2+1):N];
        end
            
        
    elseif (nargin==3)
    % truncate according to vec_trunc

        ind_trunc = vec_trunc;
    
        
    end
    
    u_hat = u_hat(ind_trunc);
       
    
end