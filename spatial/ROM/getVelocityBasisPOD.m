function [B,div_free,Vbc,rom_yM] = getVelocityBasisPOD(snapshots,sample_index,options)
% build basis for velocity field
% returns: 
% * basis B, 
% * div_free flag, 
% * non-hom BC field Vbc 
% * non-hom BC contributions yM

%% input checks
% concatenate snapshot matrices
V_total_snapshots = [snapshots.uh_total';snapshots.vh_total'];

% check input dimensions
Nspace  = size(V_total_snapshots,1); % total number of unknowns (Nu+Nv) of the original model
Nu      = options.grid.Nu;
Nv      = options.grid.Nv;
Np      = options.grid.Np;

if (Nspace ~= Nu+Nv)
    error('The dimension of the snapshot matrix does not match the input dimensions in the parameter file');
end



%% check whether snapshots are divergence free
% this gives max div for each snapshot:
div_snapshots = max(abs(options.discretization.M*V_total_snapshots + options.discretization.yM),[],1); %
% max over all snapshots:
maxdiv_snapshots = max(div_snapshots);
if (maxdiv_snapshots > 1e-14 && options.rom.rom_bc<2)
    warning(['snapshots not divergence free: ' num2str(maxdiv_snapshots)]);
end


%% subtract non-homogeneous BC contribution:

% note uh_total is stored as a Nt*Nu matrix, instead of the Nu*Nt matrix
% which we use for the SVD
Om     = options.grid.Om;
Om_inv = options.grid.Om_inv;

if (options.rom.rom_bc == 1)
    % check if the Vbc field has been stored as part of the FOM
    if (isfield(snapshots,'Vbc'))
        Vbc = snapshots.Vbc;
    else
        disp('computing Vbc field...');
        f       = options.discretization.yM;
        dp      = pressure_poisson(f,options.time.t_start,options);
        Vbc     = - Om_inv.*(options.discretization.G*dp);
    end
    V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
    %                 options.rom.yM = zeros(Np,1);
    rom_yM = zeros(Np,1);
elseif (options.rom.rom_bc == 0)
    % for rom_bc=0, we have homogeneous BC
    Vbc = zeros(Nu+Nv,1);
    %                 options.rom.yM = zeros(Np,1);
    rom_yM = zeros(Np,1);
elseif (options.rom.rom_bc == 2)
    % for rom_bc=2, we have time-dep BC and an alternative
    % method is used that uses a pressure basis
    Vbc = zeros(Nu+Nv,1);
    % store yM = -M*V
    rom_yM = -options.discretization.M*V_total_snapshots;
    if (options.rom.Mp > rank(rom_yM))
        warning('Number of pressure modes larger than rank of divergence of snapshot matrix');
    end
end



% select snapshots
V_svd = V_total_snapshots(:,sample_index);

% clear V_total_snapshots;


%% get Vbc into options (this has to be outside the j==1 if statement)
% note that the options structure gets overwritten for parametric
% studies, but Vbc and rom_yM are not overwritten
% options.rom.Vbc = Vbc;

% options.rom.yM  = rom_yM;


%% construct basis through SVD or eigenvalue problem
svd_start = toc;

% enforce momentum conservation (works for periodic domains)
if (options.rom.mom_cons == 1 && options.rom.weighted_norm == 0)
    
    e_u = zeros(Nspace,1);
    e_v = zeros(Nspace,1);
    e_u(1:Nu)     = 1;
    e_v(Nu+1:end) = 1;
    e = [e_u e_v];
    e = e / norm(e);
    
    % 1) construct (I-ee')*V_svd
    Vmod = V_svd - e*(e'*V_svd);
    % 2) take SVD
    [W,S,Z] = svd(Vmod,'econ');
    % 3) add e
    W = [e W];
    
    %     disp('error in representing vector y before truncating:');
    %     norm(Wext*Wext'*e - e,'inf')
    
elseif (options.rom.mom_cons == 1 && options.rom.weighted_norm == 1)
    
    Om_mat     = spdiags(Om,0,Nu+Nv,Nu+Nv);
    Om_sqrt    = spdiags(sqrt(Om),0,Nu+Nv,Nu+Nv);
    Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
    
    e_u = zeros(Nspace,1);
    e_v = zeros(Nspace,1);
    e_u(1:Nu)     = 1;
    e_v(Nu+1:end) = 1;
    e = [e_u e_v];
    % scale e such that e'*Om*e = I
    e = e / sqrt(norm(e'*(Om_mat*e)));
    
    % 1) construct (I-ee')*Om*V_svd
    Vmod = V_svd - e*(e'*(Om_mat*V_svd));
    % 2) apply weighting
    Vmod = Om_sqrt*Vmod;
    % 3) perform SVD
    [W,S,Z] = svd(Vmod,'econ');
    % 4) transform back
    W = Om_invsqrt*W;
    % 5) add e
    W = [e W];
    
elseif (options.rom.mom_cons == 0 && options.rom.weighted_norm == 0)
    
    % perform SVD
    %     [W,S,Z] = svd(V_svd,'econ');
    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [W,S] = getBasis(V_svd,options);
    
elseif (options.rom.mom_cons == 0 && options.rom.weighted_norm == 1)
    
    Om_sqrt    = spdiags(sqrt(Om),0,Nu+Nv,Nu+Nv);
    Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
    
    % make weighted snapshot matrix
    Vmod = Om_sqrt*V_svd;
    % perform SVD
    %     [W,S,Z] = svd(Vmod,'econ');
    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [W,S] = getBasis(Vmod,options);
    
    % transform back
    W = Om_invsqrt*W;
    
else
    error('wrong option for weighted norm or momentum conservation');
    
end

% keep V_svd in memory for case of parametric studies
% could also decide to do basis computation only for j==1
% clear V_svd;

svd_end = toc-svd_start

% take first M columns of W as a reduced basis
% maximum:
% M = size(Wu,2);
% reduction:
% M = floor(Nspace/100);
% M = 16;
M = options.rom.M;
% (better is to look at the decay of the singular values in S)
B  = W(:,1:M);
% Bu = Wu(:,1:M);
% Bv = Wv(:,1:M);
% Bu = B(1:Nu,:);
% Bv = B(Nu+1:end,:);
% options.rom.B = B;
% options.rom.Bu = Bu;
% options.rom.Bv = Bv;
% options.rom.BuT = BuT;
% options.rom.BvT = BvT;
% toc

%% relative information content:
if (size(S,2)>1)
    Sigma = diag(S);
else
    Sigma = S;
end
RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
disp(['relative energy captured by SVD = ' num2str(RIC)]);
figure(21)
semilogy(Sigma/Sigma(1),'s');
% or alternatively
% semilogy(Sigma.^2/sum(Sigma.^2),'s');


%% check whether basis is divergence free
% this gives max div for each basis vector (columns of B):
% note that yM should NOT be included here, it was already used in
% subtracting Vbc from the snapshots matrix
div_basis = max(abs(options.discretization.M*B),[],1); %
% max over all columns:
maxdiv_basis = max(div_basis);
if (options.rom.rom_bc < 2)
    if (maxdiv_basis > 1e-12)
        warning(['ROM basis not divergence free: ' num2str(maxdiv_basis) '\n']);
        answer = questdlg('Do you want to add a basis for the pressure?','Basis not divergence free','Yes','No','No');
        switch answer
            case 'Yes'
                div_free = 0;
            case 'No'
                % continue as if basis is divergence-free
                div_free = 1;
        end
    else
        div_free = 1;
    end
elseif (options.rom.rom_bc == 2)
    % unsteady boundary conditions (of normal velocity component):
    
    div_free = 0;
end

