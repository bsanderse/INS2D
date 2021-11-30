function [Breal,Breal_inv] = DFT_matrix_real(N,Ntrunc,vec_trunc,method)
% DFT_MATRIX_REAL: return 1D discrete Fourier transform basis for real data
% vectors, transforming data vector into a set of (N/2+1) real and (N/2-1) imaginary
% components
% note that the definition here is such that
% B*B_inv = B_inv*B = I, if no truncation is used

% the truncation parameter Ntrunc is used to select a subset of rows from B (cols
% of Binv), obtained by ordering in terms of frequencies.
% this is done in such a way that the resulting matrix B has roughly size Ntrunc*N
% note that each frequency has a positive and a negative component, but
% only the positive ones are selected for the real transform
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
% of the signal, or to skip positive or negative frequencies.

if (mod(N,2)~=0)
    error('input N to DFT_matrix_real should be an even integer');
end

if (nargin<4)
    % options: 'direct', 'transform_complex'
    % default is 'direct'
    method = 'direct';
end

switch method
    
    case 'direct'
        % build the DFT matrix directly
        
        %
        j  = 0:N-1; % spatial index is from 0..N-1
        k1 = 0:N/2; % frequency array for real part: constant mode, positive frequencies, Nyquist frequency
        k2 = 1:N/2-1; % frequency array for im part: positive frequencies
        
        [K1,J1] = ndgrid(k1,j);
        [K2,J2] = ndgrid(k2,j);
        
        % involve sqrt(2/N) to make Binv = B'
        Breal     = sqrt(2/N)*[cos(2*pi*(K1.*J1)/N); -sin(2*pi*(K2.*J2)/N)];

        % the resulting Breal has the following row structure:
        % index 1: constant mode (k=0)
        % index 2:N/2: cosine terms (positive freqs) (k=1:N/2)
        % index N/2+1: Nyquist (k=N/2)
        % index N/2+2:N: sine terms (positive freqs) (k=N/2+1:N-1)        
        
        % note! adapt k=0 and k=N/2:
        Breal(1,:)     = Breal(1,:)/sqrt(2);
        Breal(N/2+1,:) = Breal(N/2+1,:)/sqrt(2);
        

        
        % full expression:
        % Binv  = sqrt(2/N)*[cos(2*pi*(K1.*J1).'/N) -sin(2*pi*(K2.*J2).'/N)];
        % faster alternative:
        Breal_inv  = Breal'; % note that B is real, so ' and .' are the same
        
        % note! adapt k=0 and k=N/2:
        % Binv(:,1)     = Binv(:,1)/sqrt(2);
        % Binv(:,N/2+1) = Binv(:,N/2+1)/sqrt(2);
        
        
        
    case 'transform_complex'
        
        
        % alternative, based on transforming the complex exponential form
        
        Bcomplex = DFT_matrix(N);
        
        
        %% 1: standard approach: (not orthonormal)
        % projection from limited real coefficients to all real coefficients
        %         e1     = ones(N/2+1,1);
        %         Pr    = [spdiags(e1,0,N/2+1,N/2+1); flip(spdiags(e1,1,N/2-1,N/2+1))];
        %         Pi    = [spdiags(e1,-1,N/2+1,N/2-1); flip(spdiags(-e1,0,N/2-1,N/2-1))];
        
        % Qr is not uniquely defined
        % we take the same structure as P':
        %         e2    = e1;
        %         e2(2:end-1) = 0.5;
        %         Qr    = [spdiags(e2,0,N/2+1,N/2+1) flip(spdiags(e2(2:end),-1,N/2+1,N/2-1))];
        %         Qi    = [spdiags(e2(2:end),1,N/2-1,N/2+1) flip(spdiags(-e2(2:end),0,N/2-1,N/2-1))];
        
        % this works fine, and gives the following set of transformation matrices,
        % with the only issue that it is not orthogonal
        % Breal     = [Qr*real(Btest); Qi*imag(Btest)];
        % Breal_inv = [real(Btest)'*Pr imag(Btest)'*Pi];
        
        
        %% 2: scaled approach: (orthonormal)
        
        %         scaling    = sqrt(2)*ones(N,1);
        %         scaling(1) = 1;
        %         scaling(N/2+1) = 1;
        %         Om         = spdiags(scaling,0,N,N);
        
        % with this scaling, we get the following relations:
        % Omi       = Om(N/2+2:end,N/2+2:end);
        % Omr       = Om(1:N/2+1,1:N/2+1);
        % (Omi*Qi)' = Pi/Omi;
        % (Omr*Qr)' = Pr/Omr;
        % this means we can define
        % Qr_orth = Omr*Qr
        % Pr_orth = Pr/Omr
        
        % in this case the Q and P matrices are related
        scaling    = 0.5*sqrt(2)*ones(N,1);
        scaling(1) = 1;
        scaling(N/2+1) = 1;
        Om         = spdiags(scaling,0,N,N);
        
        % with this scaling, we get the following relations:
        Omr       = Om(1:N/2+1,1:N/2+1);
        Omi       = Om(N/2+2:end,N/2+2:end);
        e1        = ones(N/2+1,1);
        Qr_orth   = Omr*[spdiags(e1,0,N/2+1,N/2+1) flip(spdiags(e1,-1,N/2+1,N/2-1))];
        Qi_orth   = Omi*[spdiags(e1,1,N/2-1,N/2+1) flip(spdiags(-e1,0,N/2-1,N/2-1))];
        
        Pr_orth = Qr_orth';
        Pi_orth = Qi_orth';
        
        % the resulting matrices are orthonormal (Breal' = Breal_inv):
        Breal     = [Qr_orth*real(Bcomplex); Qi_orth*imag(Bcomplex)];
        Breal_inv = [real(Bcomplex)'*Pr_orth imag(Bcomplex)'*Pi_orth];
        
    otherwise
        
        error('invalid option for method in DFT_matrix_real');
end



% % truncate to Ntrunc/2 modes
% if (nargin>1 && Ntrunc>0)
%     Ntrunc    = Ntrunc/2;
% %     B         = [B(1:Ntrunc+1,:); B(end-Ntrunc+1:end,:)]; % remove rows
% %     Binv      = [Binv(:,1:Ntrunc+1) Binv(:,end-Ntrunc+1:end)];  % remove columns
%     Breal         = [Breal(1:Ntrunc+1,:); Breal(N/2+2:N/2+Ntrunc,:)]; % remove rows
%     Breal_inv      = [Breal_inv(:,1:Ntrunc+1) Breal_inv(:,N/2+2:N/2+Ntrunc)];  % remove columns
%
% end

% truncate if necessary
if (nargin>1 && (~isempty(Ntrunc) || ~isempty(vec_trunc)) )
    
    if (~isempty(Ntrunc) && Ntrunc>=N)
        warning('truncation number should be smaller or equal to matrix dimension, not truncating');
        return;
    end
    
    if (~exist('vec_trunc') ||  isempty(vec_trunc))
        % if vector is specified, we use Ntrunc for truncation
        % truncate to Ntrunc modes        
        switch method
            
            case 'direct'
                
                % index 1: constant mode (k=0)
                % index 2:N/2: cosine terms (positive freqs) (k=1:N/2)
                % index N/2+1: Nyquist (k=N/2)
                % index N/2+2:N: sine terms (positive freqs) (k=N/2+1:N-1)
                ind_trunc = [1:Ntrunc/2+1 N/2+2:N/2+Ntrunc/2+1];
                Breal     = Breal(ind_trunc,:);
                Breal_inv = Breal';
                
            case 'transform_complex'
                
                if (Ntrunc>0 && Ntrunc<N)
                    % check whether Ntrunc is odd, but leave the case of Ntrunc=1 intact
                    if (Ntrunc~=1 && mod(Ntrunc,2) ~= 0)
                        Ntrunc = 2*ceil(Ntrunc/2);
                        warning(['truncation parameter is not even, changing to Ntrunc=' Ntrunc]);
                    end
                    Bcomplex_trunc = DFT_matrix(N,Ntrunc);

                    % similar to the non-truncated case, except that the
                    % Nyquist frequency is not included and so
                    % scaling(N/2+1) does not need to be reset
                    scaling    = 0.5*sqrt(2)*ones(Ntrunc+1,1);
                    scaling(1) = 1;

                    Om         = spdiags(scaling,0,Ntrunc+1,Ntrunc+1);
                    
                    % with this scaling, we get the following relations:
                    Omi       = Om(Ntrunc/2+2:end,Ntrunc/2+2:end);
                    Omr       = Om(1:Ntrunc/2+1,1:Ntrunc/2+1);
                    % (Omi*Qi)' = Pi/Omi;
                    % (Omr*Qr)' = Pr/Omr;
                    % this means we can define
                    % Qr_orth = Omr*Qr
                    % Pr_orth = Pr/Omr
                    
                    % in this case the Q and P matrices are related
                    e1         = ones(Ntrunc/2+1,1); 
                    Qr_orth    = Omr*[spdiags(e1,0,Ntrunc/2+1,Ntrunc/2+1) flip(spdiags(e1,0,Ntrunc/2+1,Ntrunc/2))];
                    Qi_orth    = Omi*[spdiags(e1,1,Ntrunc/2,Ntrunc/2+1) flip(spdiags(-e1,0,Ntrunc/2,Ntrunc/2))];
                    
                    Pr_orth = Qr_orth';
                    Pi_orth = Qi_orth';
                    
                    % the resulting matrices are orthonormal (Breal' = Breal_inv):
                    Breal     = [Qr_orth*real(Bcomplex_trunc); Qi_orth*imag(Bcomplex_trunc)];
                    Breal_inv = [real(Bcomplex_trunc)'*Pr_orth imag(Bcomplex_trunc)'*Pi_orth];
                    % or: Breal_inv = Breal.';
                end
                
        end
        
    else
        % truncate according to vec_trunc
        %        
        Breal     = Breal(vec_trunc,:);
        Breal_inv = Breal';
                       
        
    end
    
    
end

end

