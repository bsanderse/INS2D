function u = IRDFT(u_hat_real,N,vec_trunc)
% IRDFT: return 1D real inverse discrete Fourier transform
% the definition here is such that the transform is orthonormal

% the 'extension' parameter N is used to map from a low number of modes
% N_hat to a larger number of spatial points N

% the extension vector vec_trunc can be used instead of N to indicate
% which rows were selected when computing u_hat_real, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies

N_hat     = size(u_hat_real,1);
N_vec = size(u_hat_real,2);
I = sqrt(-1);


% extend the transformed value u_hat
if (nargin>1)
    
    u_hat = zeros(N,N_vec);    

    
    if (nargin==2)
    % truncate to N_hat modes
        
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
%             ind_trunc = [1:N/2+1 (N-N/2+1):N];

            % in this case N_hat is typically odd, and Nyquist frequency is missing
            ind_cos  = 1:floor(N_hat/2)+1;
            ind_sin  = floor(N_hat/2)+2:N_hat;
            
            % get indices in the extended vector
            ind_pos = ind_cos;
            ind_neg = N - floor(N_hat/2) + 1:N;
        end
        
        
    elseif (nargin==3)
    % extend according to vec_trunc
        
        ind_cos  = 1:floor(N_hat/2)+1;
        ind_sin  = floor(N_hat/2)+2:N_hat;
        
        ind_pos = vec_trunc(vec_trunc<=N/2);
        ind_neg =  flip(N-vec_trunc(vec_trunc<=N/2 & vec_trunc>1)+2);% vec_trunc(vec_trunc<=N/2 & vec_trunc>1);
        
        
    end
          
    % scale back the sqrt(2) term that was introduced in the forward
    % transform; note: not for k=0 
    u_hat_real(ind_cos(2:end),:) = u_hat_real(ind_cos(2:end),:)/sqrt(2);
    u_hat_real(ind_sin,:) = u_hat_real(ind_sin,:)/sqrt(2); 


    % set real component, positive frequencies
    u_hat(ind_pos,:) = u_hat_real(ind_cos,:);
    % add imaginary component, positive frequencies
    u_hat(ind_pos(2:end),:) = u_hat(ind_pos(2:end),:) + I*u_hat_real(ind_sin,:);

    % for negative frequencies, we have to 'flip'
    % set real component, negative frequencies
    u_hat(ind_neg,:) = u_hat_real(flip(ind_cos(2:end)),:);
    % add imaginary component, negative frequencies
    u_hat(ind_neg,:) = u_hat(ind_neg,:) - I*u_hat_real(flip(ind_sin),:);

    % then do the ifft on the complex vector
    u = ifft(u_hat)*sqrt(N);    
    
else
    
    % transform back to complex values
    ind_cos  = 1:floor(N_hat/2)+1;
    ind_sin  = floor(N_hat/2)+2:N_hat;
    ind_pos  = ind_cos;
    ind_neg  = ind_sin;
    u_hat = zeros(N_hat,N_vec);
    
    % scale back the sqrt(2) term that was introduced in the forward
    % transform; note: not for k=0 and k=N/2
    u_hat_real(ind_cos(2:end-1),:) = u_hat_real(ind_cos(2:end-1),:)/sqrt(2);
    u_hat_real(ind_sin,:) = u_hat_real(ind_sin,:)/sqrt(2);
   
    % set real component, positive frequencies
    u_hat(ind_pos,:) = u_hat_real(ind_cos,:);
    % add imaginary component, positive frequencies
    u_hat(ind_pos(2:end-1),:) = u_hat(ind_pos(2:end-1),:) + I*u_hat_real(ind_sin,:);

    % for negative frequencies, we have to 'flip'
    % set real component, negative frequencies
    u_hat(ind_neg,:) = u_hat_real(flip(ind_cos(2:end-1)),:);
    % add imaginary component, negative frequencies
    u_hat(ind_neg,:) = u_hat(ind_neg,:) - I*u_hat_real(flip(ind_sin),:);

    % then do the ifft on the complex vector
    u = ifft(u_hat)*sqrt(N_hat);
    
end

end