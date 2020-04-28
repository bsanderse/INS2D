%% test method of snapshots (MOS), and compare with SVD and direct method

clear all
close all

%%
% select which algorithms to test
% note: testing the direct method can be extremely expensive for large problems
svd_test = 1;
direct_test = 0;
ss_test = 1;

snapshot_data = 'results/LDC_unsteady_1.000e+02_20x20_FE/matlab_data.mat';
% snapshot_data = 'results/actuator_ROM_unsteady_force/matlab_data.mat';
% snapshot_data = 'results/LDC_unsteady_rerun_April2020/matlab_data.mat';
% snapshot_data = 'results/shear_layer_ROM_snapshots_rerunApril2020/matlab_data.mat';


snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','options');
%
X = [snapshots.uh_total';snapshots.vh_total'];
% X = X - mean(X); %snapshots.options.rom.V_bc;
[Nv,Ns] = size(X);

% number of modes:
M = snapshots.options.rom.M;

%%
% should we subtract the mean?
% V_svd = V_total_snapshots - mean(V_total_snapshots);
% V_svd = rand(20000,200);

%% SVD
if (svd_test == 1)
    % standard SVD on snapshot matrix, which is NV*Nsnapshots
    % this is preferred for NV<Nsnapshots
    % note that W are the eigenvectors of V_svd*V_svd'
    tic
    svd_start   = toc;
    % for Nv>>Ns, we get for the economic SVD:
    % Phi_svd: Nv*Ns
    % S_svd: Ns*Ns
    % Psi_svd: Ns*Ns
    [Phi_svd,S_svd,Psi_svd] = svd(X,'econ');
    S_svd       = diag(S_svd); % singular values
    d_svd       = S_svd.^2; % eigenvalues
    svd_time    = toc - svd_start
    
    
    figure
    semilogy(S_svd,'s-');
end

%% direct method
if (direct_test==1)
    direct_start = toc;
    R     = X*X';
    % eigenvectors correspond to basis Phi
    % Phi: Nv*Nv
    % [Phi_direct,d_direct] = eig(R,'vector');
    [Phi_direct,D_direct] = eigs(R,M,'largestreal');
    d_direct = diag(D_direct);
    % sort the eigenvectors according to descending eigenvalues
    [d_sort,ind]    = sort(d_direct,'descend');
    Phi_direct_sort = Phi_direct(:,ind);
    direct_time = toc - direct_start
    
end

%% method of snapshots:
% see e.g.
% https://depositonce.tu-berlin.de/bitstream/11303/9456/5/podnotes_aiaa2019.pdf
% https://web.stanford.edu/group/frg/course_work/CME345/CA-CME345-Ch4.pdf
% http://people.math.sc.edu/wangzhu/reference/Wang2016approximate.pdf

if (ss_test == 1)
    % get eigenvectors of V_svd'*V_svd
    snapshot_start = toc;
    C     = X'*X;
    % kappa = cond(C);
    % r     = rank(C);
    
    % eigenvectors correspond to basis Psi
    % Psi: Ns*Ns
    % full eig routine:
%     [Psi_ss,d_snapshots] = eig(C,'vector');
%     [d_sort,ind] = sort(d_snapshots,'descend');
    % sort the eigenvectors according to descending eigenvalues
%     Psi_ss_sort = Psi_ss(:,ind);
    % remove zero eigenvalues and truncate
%     r = min(rank(C),M);
%     d_sort      = d_sort(1:r,1);
%     Psi_ss_sort = Psi_ss_sort(:,1:r);


    % eigs routine, already orders the eigenvalues:
    [Psi_ss_sort,D_snapshots] = eigs(C,M,'largestreal');
    d_ss = diag(D_snapshots);
            
    %
    % r       = length(d_sort);
    
    % get Phi based on Psi
    % this requires that Psi'*Psi = I_r
    disp('orthogonality of right singular vectors of snapshot matrix:')
    norm(Psi_ss_sort'*Psi_ss_sort - speye(M),Inf)
    % note: the division through the square root of small eigenvalues can cause
    % issues in numerical accuracy
    S_ss        = spdiags(1./sqrt(d_ss),0,M,M);
    Phi_ss_sort = (X*Psi_ss_sort)*S_ss;
    
    snapshot_time = toc-snapshot_start
    
end

%% compare approaches
r = M;
% r = min(rank(C),rank(R));

if (svd_test == 1)
    W_svd    = Phi_svd(:,1:r);
end
if (direct_test == 1)
    W_direct = Phi_direct_sort(:,1:r);
end
if (ss_test == 1)
    W_ss     = Phi_ss_sort(:,1:r);
end

% normalize column of each matrix
% W4 = W2;
% W4 = normc(W2);

eps=1e-14;
% NOTE: the modes calculated separately by the two methods may have opposite
% signs
norm_svd_direct = zeros(r,1);
norm_ss_direct  = zeros(r,1);
norm_svd_ss     = zeros(r,1);

for i=1:r-1
    if (svd_test == 1 && direct_test == 1)
        norm1 = norm(W_svd(:,i)-W_direct(:,i),Inf);
        norm2 = norm(W_svd(:,i)+W_direct(:,i),Inf);
        norm_svd_direct(i) = min(norm1,norm2);
    end

    if (ss_test == 1 && direct_test == 1)
        
        norm1 = norm(W_ss(:,i)-W_direct(:,i),Inf);
        norm2 = norm(W_ss(:,i)+W_direct(:,i),Inf);
        norm_ss_direct(i) = min(norm1,norm2);
    end

    if (svd_test == 1 && ss_test == 1)
        
        norm1 = norm(W_svd(:,i)-W_ss(:,i),Inf);
        norm2 = norm(W_svd(:,i)+W_ss(:,i),Inf);
        norm_svd_ss(i) = min(norm1,norm2);
    end
    
end

% semilogy(abs(sqrt(d_sort)-S_diag(1:r)))
figure
semilogy(norm_svd_direct,'s-')
hold on
semilogy(norm_ss_direct,'o-')
semilogy(norm_svd_ss,'d-')
legend('SVD vs direct','Snapshots vs direct','SVD vs snapshots');
ylabel('error between different algorithms')


% check orthogonality:
if (svd_test == 1)
    disp('orthogonality of SVD basis:')
    norm((W_svd'*W_svd)-speye(r),Inf)
end
if (direct_test == 1)
    disp('orthogonality of direct basis:')
    norm((W_direct'*W_direct)-speye(r),Inf)
end
if (ss_test == 1)
    disp('orthogonality of snapshot basis:')
    norm((W_ss'*W_ss)-speye(r),Inf)
    % if this value is not close to zero, but the value of     
    % norm(Psi_ss_sort'*Psi_ss_sort - speye(M),Inf) IS close to zero,
    % then the issue lies in the accuracy of the computation of 1./eigenvalues

end

%% check divergence freeness of the different approaches and projection error of the basis
Div = snapshots.options.discretization.M; % size Np*Nv
maxdiv_svd = zeros(r,1);
maxdiv_direct = zeros(r,1);
maxdiv_ss = zeros(r,1);
proj_error_svd = zeros(r,1);
proj_error_direct = zeros(r,1);
proj_error_ss = zeros(r,1);
R = X*X';
for i=1:r-1
    
    if (svd_test == 1)
        maxdiv_svd(i)    = max(abs(Div*W_svd(:,i)));
        proj_error_svd(i) = norm( R*W_svd(:,i) - d_svd(i)*W_svd(:,i));
    end
    if (direct_test == 1)
        maxdiv_direct(i) = max(abs(Div*W_direct(:,i)));
        proj_error_direct(i) = norm( R*W_direct(:,i) - d_direct(i)*W_direct(:,i));
    end
    if (ss_test == 1)
        maxdiv_ss(i)     = max(abs(Div*W_ss(:,i)));
        proj_error_ss(i) = norm( R*W_ss(:,i) - d_ss(i)*W_ss(:,i));

    end
end
figure
semilogy(maxdiv_svd,'s-');
hold on
semilogy(maxdiv_direct,'o-');
semilogy(maxdiv_ss,'d-');
legend('SVD','Direct','Snapshots');
ylabel('max divergence');

figure
semilogy(proj_error_svd,'s-');
hold on
semilogy(proj_error_direct,'o-');
semilogy(proj_error_ss,'d-');
legend('SVD','Direct','Snapshots');
ylabel('projection error');

% any errors should come from the fact that X*X^T Phi does not equal Phi
% accurately enough:
