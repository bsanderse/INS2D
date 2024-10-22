% todo: make function out of this!

switch options.rom.rom_type
    
    case 'POD'
        
        [V_svd,Vbc,snapshots] = load_snapshot_data(snapshot_data,options,dt_sample,t_sample);
        
        %% get Vbc into options (this has to be outside the j==1 if statement)
        options.rom.Vbc = Vbc;
        
        
        %% construct basis through SVD or eigenvalue problem
        svd_start = toc;

        Om = options.grid.Om;
        Nu      = options.grid.Nu;
        Nv      = options.grid.Nv;
        
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

        %% store operator inference snapshot data
        if options.rom.opinf_type ~= "intrusive" &&  ~exist("opinf_snapshot_data",'var')

            opinf_V_svd = V_svd;

            [A,A_dot] = get_opinf_snapshots(B'*(Om.*opinf_V_svd),snapshots.dt); % presumably only works for rom_type == "POD"
            
            %% extend snapshot data by 90 degree-rotated data (modulo minus sign)
            % Nu = options.grid.Nu;
            % opinf_V_svd_rot = [V_svd(Nu+1:end,:); V_svd(1:Nu,:)];
            % [A_rot,A_dot_rot] = get_opinf_snapshots(B'*(Om.*opinf_V_svd_rot),snapshots.dt); % presumably only works for rom_type == "POD"
            % A = [A A_rot];
            % A_dot = [A_dot A_dot_rot];
            %%
            
            options.rom.A = A;
            options.rom.A_dot = A_dot;
        end
        %%

        clear V_svd;
        
        svd_end(j) = toc-svd_start
       
        
        % relative information content:
        if (size(S,2)>1)
            Sigma = diag(S);
        else
            Sigma = S;
        end
        RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
        disp(['relative energy captured by SVD = ' num2str(RIC)]);
        figure
        semilogy(Sigma/Sigma(1),'s');
        % or alternatively
        % semilogy(Sigma.^2/sum(Sigma.^2),'s');
        
        %% check whether basis is divergence free
        % this gives max div for each basis vector (columns of B):
        % note that yM should not be included here, it was already used in
        % subtracting Vbc from the snapshots matrix
        div_basis = max(abs(options.discretization.M*B),[],1); %
        % max over all columns:
        maxdiv_basis = max(div_basis);
        if (maxdiv_basis > 1e-14)
            warning(['ROM basis not divergence free: ' num2str(maxdiv_basis) '\n']);
            % warning('Adding basis for pressure\n');
            % options.rom.div_free = 0;
            options.rom.div_free = 1;
        else     
            options.rom.div_free = 1;
        end
        
        %% pressure recovery (postprocessing)
        if (options.rom.pressure_recovery == 1)
            disp('computing SVD of pressure snapshots...');
            svd_start2 = toc;
            % note p_total is stored as a Nt*Np matrix instead of Np*Nt which we use for
            % the SVD
            % use same snapshot_indx that was determined for velocity
            
            % select snapshots
            p_total_snapshots = snapshots.p_total';
            if (options.rom.pressure_mean == 1)
                % subtract temporal mean
                options.rom.p_mean = mean(p_total_snapshots,2);
                p_total_snapshots = p_total_snapshots - options.rom.p_mean;
            end
            p_svd  = p_total_snapshots(:,snapshot_sample);
            
            % take first Mp columns of Wp as a reduced basis
            % (better is to look at the decay of the singular values in Sp to determine M)
            if (isfield(options.rom,'Mp'))
                Mp = options.rom.Mp;
            else
                % if not defined, use same number of modes as for velocity
                warning('number of pressure modes not defined, defaulting to number of velocity modes');
                Mp = options.rom.M;
            end
            
            if (options.rom.weighted_norm == 0)
                
                [Wp,Sp] = getBasis(p_svd,options,options.rom.Mp);
                
                % perform SVD
                %     [Wp,Sp,Zp] = svd(p_svd,'econ');
                
            elseif (options.rom.weighted_norm == 1)
                
                Np          = options.grid.Np;
                Omp         = options.grid.Omp;
                Omp_sqrt    = spdiags(sqrt(Omp),0,Np,Np);
                Omp_invsqrt = spdiags(1./sqrt(Omp),0,Np,Np);
                
                % make weighted snapshot matrix
                pmod = Omp_sqrt*p_svd;
                
                % getBasis can use different methods to get basis: SVD/direct/snapshot
                % method
                [Wp,Sp] = getBasis(pmod,options);
                
                % transform back
                Wp = Omp_invsqrt*Wp;
            end
            
            Bp = Wp(:,1:Mp);
            options.rom.Bp = Bp;
            
            svd_end(j) = svd_end(j) + toc - svd_start2
            
            hold on
            if (size(Sp,2)>1)
                SigmaP = diag(Sp);
            else
                SigmaP = Sp;
            end
            semilogy(SigmaP/SigmaP(1),'o');
            
        end
        
        
        
    case 'Fourier'
        
        Nu  = options.grid.Nu;
        Nv  = options.grid.Nv;

        NV  = options.grid.NV;
        Vbc = zeros(NV,1);
        options.rom.Vbc = Vbc;

        % construct velocity basis B from discrete Fourier transform
        % note that we have V = B*R, which we split into two components
        % u = B_u*R_u, v = B_v*Rv
        Phi = getFourierBasis(options); %,M);
        % scale with sqrt(Omega) in order to get Phi'*Om*Phi = I
        Phi = 1/sqrt(options.grid.Om(1)) * Phi;
        B   = [Phi spalloc(Nu,M,0); ...
               spalloc(Nv,M,0) Phi];
        options.rom.B = B;
        M   = size(B,2);
        options.rom.M = M;

        % we also need to make a pressure basis Bp because B is not
        % divergence-free
%         Bp = B;
        options.rom.Bp = Phi;
        
        % set pressure_recovery to 2 in order that reduced Poisson operator
        % is constructed
        div_basis = max(abs(options.discretization.M*B),[],1); %
        % max over all columns:
        maxdiv_basis = max(div_basis);
        if (maxdiv_basis > 1e-14)
            disp(['Fourier basis is not divergence free: ' num2str(maxdiv_basis)]);
            disp('Adding basis for pressure');
            options.rom.div_free = 0;
        else     
            options.rom.div_free = 1;
        end
        
%        
%         % make the basis divergence free
%         Om_inv = options.grid.Om_inv;
%         G  = options.discretization.G;
%         for i=1:M
%            
%             % make velocity field divergence free            
%         
%             f  = options.discretization.M*B(:,i) + options.discretization.yM;
%             dp = pressure_poisson(f,t,options);
%             
%             % look at the resulting modes, before and after projection
% %             figure(1)
% %             surf(reshape(B(options.grid.indu,i),Nx,Ny)')
%             
%             B(:,i)  = B(:,i) - Om_inv.*(G*dp);
%             % turns out that they effectively become 0!
% %             figure(2)
% %             surf(reshape(B(options.grid.indu,i),Nx,Ny)')
%         end
%         
%         % make orthogonal
%         [Bnew,R] = qr(B);
% %         B = Bnew/sqrt(options.grid.Om(1));
%         
%         [B,R] = mwgson(B,options.grid.Om);
        
        div_basis = max(abs(options.discretization.M*B),[],1); %
        % max over all columns:
        maxdiv_basis = max(div_basis);
        if (maxdiv_basis > 1e-14)
            disp(['Fourier basis is not divergence free: ' num2str(maxdiv_basis)]);
            disp('Adding basis for pressure');
            options.rom.div_free = 0;
        else     
            options.rom.div_free = 1;
        end
        
        disp('orthogonality of reduced basis B wrt omega:');
        Om     = options.grid.Om;
        Om_mat = spdiags(Om,0,Nu+Nv,Nu+Nv);
        B_orth = max(max(abs(B'*Om_mat*B-speye(size(B,2)))))
        if (B_orth>1e-13)
            warning('Reduced basis B not orthogonal wrt Omega');
        end
        
        keyboard
        
    case 'FDG'
        
        % Finite Difference Galerkin = reformulation in streamfunction form
        
%         % generate null vectors
%         Null_u = zeros(options.grid.Nu,1);
%         Null_u(2)    = 1;
%         Null_u(2+Nx) = -1;
%         Null_v = zeros(options.grid.Nv,1);
%         Null_v(1+Nx) = -1;
%         Null_v(2+Nx) = 1;

        % use vorticity operator to construct null space
        % note: need to get rid of extra periodic entries
        
        %% for periodic BC, not covering entire mesh
        
        Nu  = options.grid.Nu;
        Nv  = options.grid.Nv;
        
        % du/dy, like Su_uy
        % note that we use hy here, and not gyd; this is to make sure that
        % W' is a null vector of M
        diag1  = 1./options.grid.hy;
        W1D    = spdiags([-diag1 diag1],[-1 0],Ny,Ny);
        W1D(1,end) = -W1D(1,1);
        % extend to 2D
        Wu_uy  = kron(W1D,speye(Nx));

        % dv/dx, like Sv_vx
        % note that we use hx here, and not gxd; this is to make sure that
        % W' is a null vector of M
        diag1  = 1./options.grid.hx;
        W1D    = spdiags([-diag1 diag1],[-1 0],Nx,Nx);
        W1D(1,end) = -W1D(1,1);
        % extend to 2D
        Wv_vx  = kron(speye(Ny),W1D);

        % curl operator that acts on velocity fields
        W = [-Wu_uy Wv_vx];
        % null space is then given by W'
        B = W';
        % is the divergence of B indeed zero?

        div_basis = max(abs(options.discretization.M*B),[],1); %
        % max over all columns:
        maxdiv_basis = max(div_basis);        
        if (maxdiv_basis > 1e-14)
            warning(['ROM basis not divergence free: ' num2str(maxdiv_basis) '\n']);
            warning('Adding basis for pressure\n');
            options.rom.div_free = 0;
        else     
            options.rom.div_free = 1;
        end            
        
        
        % we can thus expand V, such that M*V=M*B*psi=0 for any psi:
        % V = B*psi
        % and
        % psi = (B'*Om*B)^{-1} * B^T * Om * V
                
        Om     = options.grid.Om;
        Om_mat = spdiags(Om,0,Nu+Nv,Nu+Nv);
       
        % ROM dimension (=NV - Np)
        % there is no real dimension reduction as in the case of
        % truncation, it's just going from velocity to streamfunction space
        M = size(B,2);
      
      
        % get the oblique projection (since B is not
        % orthonormal)
        options.rom.B_inv = decomposition(B.'*Om_mat*B);
        

        options.rom.B = B;
        options.rom.M = M;

        Vbc = zeros(Nu+Nv,1);
        options.rom.Vbc = Vbc;
        
        
    case 'FDG-Fourier'
        
        %% for periodic BC, not covering entire mesh
        
        Nu  = options.grid.Nu;
        Nv  = options.grid.Nv;
        
        % du/dy, like Su_uy
        % note that we use hy here, and not gyd; this is to make sure that
        % W' is a null vector of M
        diag1  = 1./options.grid.hy;
        W1D    = spdiags([-diag1 diag1],[-1 0],Ny,Ny);
        W1D(1,end) = -W1D(1,1);
        % extend to 2D
        Wu_uy  = kron(W1D,speye(Nx));

        % dv/dx, like Sv_vx
        % note that we use hx here, and not gxd; this is to make sure that
        % W' is a null vector of M
        diag1  = 1./options.grid.hx;
        W1D    = spdiags([-diag1 diag1],[-1 0],Nx,Nx);
        W1D(1,end) = -W1D(1,1);
        % extend to 2D
        Wv_vx  = kron(speye(Ny),W1D);

        % curl operator that acts on velocity fields
        W = [-Wu_uy Wv_vx];
%         M_Null = W';
        % null space is then given by W'
        C = W';
        % is the divergence of C indeed zero?
        disp('divergence of curl matrix:')
        max(max(abs(options.discretization.M*C)))
        max(max(abs(C'*options.discretization.G)))
        % we can thus expand V, such that M*V=M*C*psi=0 for any psi
        % V = C*psi
        
        % now build a Fourier basis for psi (streamfunction)
        % psi = Phi*R,
        % where R are Fourier coefficients for streamfunction,
        % and Phi is the INVERSE Fourier transform
        Phi = getFourierBasis(options,M);
        
        M   = size(Phi,2);

        disp(['reduced size M = ' num2str(M)]);
       
        % in all cases, we get the final basis as follows
        % split the Omega contribution symmetrically
        
        Om     = options.grid.Om;
%         Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
        Om_mat     = spdiags(Om,0,Nu+Nv,Nu+Nv);
%         if (options.rom.weighted_norm==1)
%             C   = sqrt(hx(1)*hy(1))*C;
%         end
%         C   = Om_invsqrt*C;
        
        % then the final basis is given by
        % V = C*psi = C*Phi*R = B*R
        B   = C*Phi;
       
        
        % is B' * Om_mat * B diagonal?
        disp('is B^T * Om * B diagonal?');
        mass_matrix = B'*Om_mat*B;
        diag_mass_matrix = diag(mass_matrix);
        B_diag = max(max(abs(mass_matrix - diag(diag_mass_matrix))))
        if (B_diag>1e-12)
            warning('not diagonal');
        end
        % check for zero eigenvalues
        if (min(abs(diag_mass_matrix)) < 1e-12)
            warning('zero eigenvalue in diagonal matrix');
        end
        
        % an alternative construction that does not use the explicit
        % formation of B is as follows:
        % we want to compute Phi^T*C^T*Om*C*Phi = Phi^T*Lpsi*Psi
        % where Lpsi is the Poisson matrix for the streamfunction:
        Lpsi  = C'*Om_mat*C;
        % now convert to Fourier domain
        % first left-multiply with Phi^T, the forward FFT transform
        Lpsi_hat = zeros(M,size(Lpsi,2));
        for i=1:size(Lpsi,2)
             LFourier = RDFT2(reshape(full(Lpsi(:,i)),Nx,Ny).',2,2).'; 
             Lpsi_hat(:,i) = LFourier(:);
        end
        % now right-multiply, Lpsi_hat*Psi, which is rewritten as
        % (Phi^T * Lpsi_hat^T)^T = (RDFT(Lpsi_hat^T)^T)
        % this leads to the following "mass matrix"
        Mpsi_hat = zeros(M,M);
        for i=1:size(Lpsi_hat,1)
             LFourier = RDFT2(reshape(full(Lpsi_hat(i,:)),Nx,Ny).',2,2).'; 
             Mpsi_hat(:,i) = LFourier(:);
        end
       
        % with these routines, we can also build B alternatively as follows
        Balt = zeros(size(B));
        for i=1:size(C,1)
            CFourier = RDFT2(reshape(full(C(i,:)),Nx,Ny).',2,2).'; 
            Balt(i,:) = CFourier(:);
        end
        
        max(max(abs(mass_matrix-Mpsi_hat)))
        max(max(abs(B-Balt)))

        
        % get the oblique projection (for the case that B is not
        % orthonormal)
%         options.rom.B_inv = decomposition(B.'*Om_mat*B);
        
        % get nullspace of div-operator (all vectors v such that M*v=0)
%         NullSpace = null(full(options.discretization.M));
        % construct basis from linear combination of the nullspace,
        % 
%         B = NullSpace(:,2:end);

        options.rom.B = B;
        options.rom.M = M;
%         options.rom.C = C;
%         options.rom.Phi = Phi;

        Vbc = zeros(Nu+Nv,1);
        options.rom.Vbc = Vbc;  
        
        options.rom.B_inv = 1./diag_mass_matrix; %(abs(mass_matrix)>1e-8); %spdiags(mass_matrix(abs(mass_matrix)>1e-8),0,options.grid.Np,options.grid.Np-1);
        
        options.rom.div_free = 1;

        
        
    otherwise
        error ('wrong ROM type')
        
end