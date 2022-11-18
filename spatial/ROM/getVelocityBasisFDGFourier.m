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