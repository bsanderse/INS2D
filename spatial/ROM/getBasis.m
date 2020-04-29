function [Phi,S] = getBasis(X,options,M)
% GETBASIS: get ROM spatial basis
% input: snapshot matrix X, size Nv*Nsnapshots; number of modes
% output: spatial basis Phi, and singular values S

basis_type = options.rom.basis_type;
if (nargin<3)
    M = options.rom.M;
end

[Nv,Nsnapshots] = size(X);

% no basis type defined, use dimensions of S to decide on method
if (basis_type == 0)
    if (Nv>Nsnapshots)
        basis_type = 3; % snapshot method
    else
        basis_type = 1; % SVD
    end
end

% get only largest eigenvalues:
use_eigs = 1;

switch basis_type
    
    case 1 %'SVD'
        disp('ROM basis constructed using SVD...');
        
        [Phi,S,~] = svd(X,'econ');
        S         = diag(S);
        
    case 2 %'direct'
        disp('ROM basis constructed using direct method...');
        
        % eigenvectors correspond to basis Phi
        % Phi: Nv*Nv
        R     = X*X';
        
        if (use_eigs == 0)
            % approach 1: use full eig computation
            [Phi_direct,d_direct] = eig(R,'vector');
            [d_sort,ind]    = sort(d_direct,'descend');
            Phi   = Phi_direct(:,ind);
            Phi   = Phi(:,1:M);
            S     = sqrt(d_sort(1:M));
            
        elseif (use_eigs == 1)
            % approach 2:
            % R has positive real eigenvalues, as it is SPD
            % eigs routine, already orders the eigenvalues:
            % only get M largest eigenvalues:
            [Phi,D_direct] = eigs(R,M,'largestreal');
            S = diag(D_direct);
        end
        
        
    case 3 % 'snapshots'
        disp('ROM basis constructed using snapshot method...');
        
        % method of snapshots, a la Sirovich
        C     = X'*X;
        
        % eigenvectors correspond to basis Psi
        % Psi: Ns*Ns
        
        if (use_eigs == 0)
            % approach 1: use full eig computation
            [Psi_ss,d_snapshots] = eig(C,'vector');
            [d_sort,ind] = sort(d_snapshots,'descend');
            %         sort the eigenvectors according to descending eigenvalues
            Psi_ss_sort = Psi_ss(:,ind);
            % truncate
            d_sort      = d_sort(1:M);
            Psi_ss_sort = Psi_ss_sort(:,1:M);
        elseif (use_eigs == 1)
            % approach 2: use eigs
            % C has positive real eigenvalues, as it is SPD
            % eigs routine, already orders the eigenvalues:
            % only get M largest eigenvalues:
            [Psi_ss_sort,D_snapshots] = eigs(C,M,'largestreal');
            d_sort = diag(D_snapshots);
        end
        
        % warn potential issues if rank of C is lower than M
        r = rank(C);
        if (r<M)
            warning('rank of C is lower than requested number of modes'); %, reverting to SVD');
            %             options.rom.basis_type = 1;
            %             [Phi,S] = getBasis(X,options);
        end
        
        
        % get Phi based on Psi
        % this requires that Psi'*Psi = I_r
        %         norm(Psi_ss_sort'*Psi_ss_sort - speye(r),Inf)
        % note: the division through the square root of small eigenvalues can cause
        % issues in numerical accuracy
        % also note that the eigenvalues can be less accurate than the
        % singular values
        S           = sqrt(d_sort);
        S_ss        = spdiags(1./S,0,M,M);
        Phi         = (X*Psi_ss_sort)*S_ss;
        
        
        
    otherwise
        error('wrong basis type');
end

end
