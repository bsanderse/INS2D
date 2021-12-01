function u = IDFT(u_hat,N,vec_trunc)
% IDFT: return 1D inverse discrete Fourier transform 
% the definition here is such that the transform is orthonormal

% the 'extension' parameter N is used to map from a low number of modes
% N_hat to a larger number of spatial points N

% the extension vector vec_trunc can be used instead of N to indicate
% which rows are to be selected, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies

N_hat = size(u_hat,1);
N_vec = size(u_hat,2);

% extend the transformed value u_hat
if (nargin>1)
    
    if (isempty(N))
        error('N needs to specified');
    end

    u_hat_full = zeros(N,N_vec);    

    if (nargin==2)
    % extend to N modes

    
        if (N<=0)
            error(['check value for parameter N, N=' N]);
        elseif (N<=N_hat)
            warning('extension parameter smaller than or equal to N_hat, no extension');
            ind_trunc = 1:N_hat;
        else
            % check whether N is odd, but leave the case of N=1 intact
            if (N~=1 && mod(N,2) ~= 0)
                N = 2*ceil(N/2);
                warning(['extension parameter is not even, changing to N=' N]);
            end
            
            % single mode: ind = 1
            % two modes: ind = [1:2, N]
            % three modes: ind = [1:3, N-1:N]
            % typically N_hat will be odd
            ind_trunc = [1:floor(N_hat/2)+1 (N-floor(N_hat/2)+1):N];
        end
            
        
    elseif (nargin==3)
    % extend according to vec_trunc

        ind_trunc = vec_trunc;
    
        
    end
    u_hat_full(ind_trunc,:) = u_hat;
    
else
    N = N_hat;
    u_hat_full = u_hat;
end

u = ifft(u_hat_full)*sqrt(N);

end