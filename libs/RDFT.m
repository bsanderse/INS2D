function u_hat_real = RDFT(u,N_hat,vec_trunc)
% RDFT: return 1D real discrete Fourier transform
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

N     = size(u,1);
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
        % TODO: if ind_trunc does not include index 1, we should check the
        % sqrt(2) scaling
        if (isempty(find(vec_trunc==1,1)))
            warning('vec_trunc does not include index 1, check scaling');
        end
        ind_trunc = vec_trunc;
        
        
    end
    
    u_hat = u_hat(ind_trunc,:);
    
    N_hat = size(u_hat,1);
    
    % N_hat is typically odd, and Nyquist frequency will be missing
    ind_cos  = 1:floor(N_hat/2)+1;
    ind_sin  = 2:floor(N_hat/2)+1;
    
    % introduce sqrt(2) scaling for orthonormaltiy
    u_hat_cos  = real(u_hat(ind_cos,:));
    u_hat_cos(2:end,:) = u_hat_cos(2:end,:)*sqrt(2);
    u_hat_sin  = imag(u_hat(ind_sin,:))*sqrt(2);
    u_hat_real = [u_hat_cos; u_hat_sin];
    
else  % only 1 argument to function call     

    ind_cos  = 1:N/2+1;
    ind_sin  = 2:N/2;
    u_hat_cos = real(u_hat(ind_cos,:));
    
    % introduce sqrt(2) scaling for orthonormaltiy
    u_hat_cos(2:end-1,:) = u_hat_cos(2:end-1,:)*sqrt(2);
    u_hat_sin  = imag(u_hat(ind_sin,:))*sqrt(2);
    
    u_hat_real = [u_hat_cos; u_hat_sin];
    
end

end