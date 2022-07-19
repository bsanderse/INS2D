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

    W0 = W;
    
    % transform back
    W = Om_invsqrt*W;

    %% mission: phi consistent phi bc
    if options.rom.rom_bc == 2
        t_js = t_start:dt:t_end;
    else
        t_js = 0;
        Mbc = 1;
    end
    X_bc = zeros(length(get_bc_vector_yBC(0,options)),length(t_js));
    for jj=1:length(t_js)
        t_j = t_js(jj);
        X_bc(:,jj) = get_bc_vector_yBC(t_j,options);
    end

    if options.verbosity.equivalence_cheat == 2
        norm(W0(:,1:M)-Vmod*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))
        norm(W0(:,1:M)-Om_sqrt*V_svd*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))
        norm(Om_invsqrt*W0(:,1:M)-V_svd*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))
        W_bc = X_bc*Vmod'*W0(:,1:Mbc)*diag(S(1:Mbc).^-2);

        %     phi_bc = W_bc(:,1:Mbc); % superfluous
        [phi_bc,R_bc] = qr(W_bc,0);

        options.rom.phi_bc = phi_bc;
        for jj = 1:Mbc
            yBC = phi_bc(:,jj);
            Y_M(:,jj) = get_yM(options,yBC);
        end
    end

    %         cond_fac = 10^-6;
%         X_bc_rank = sum(abs(Sigma_bc/Sigma_bc(1))>cond_fac);
%         if Mbc > X_bc_rank
%             Mbc = X_bc_rank;
%         end
   %%
    
else
    error('wrong option for weighted norm or momentum conservation');
    
end
% clear V_svd;

svd_end(j) = toc-svd_start

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
Bu = B(1:Nu,:);
Bv = B(Nu+1:end,:);
options.rom.B = B;
options.rom.Bu = Bu;
options.rom.Bv = Bv;
% options.rom.BuT = BuT;
% options.rom.BvT = BvT;
toc

% relative information content:
if (size(S,2)>1)
    Sigma = diag(S);
else
    Sigma = S;
end
% RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
% disp(['relative energy captured by SVD = ' num2str(RIC)]);
disp('commented some information out')

if (options.visualization.show_sigmas == 1)
    figure(123)
    semilogy(Sigma/Sigma(1),'s','displayname', 'singular values velocity snapshot matrix');
    hold on
end

% or alternatively
% semilogy(Sigma.^2/sum(Sigma.^2),'s');