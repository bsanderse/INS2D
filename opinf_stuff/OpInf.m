function [Diff, Conv] = OpInf(A_dot,A_hat,rom_type) 

    %% botch: SVD compression
        if true
            N_hat = size(A_hat,1);
            [U,S,V] = svd(A_hat',"econ");

            rank_ = sum(abs(diag(S))>sqrt(eps))
            rank_ = min(rank_,N_hat); % N tilde does not need to be larger than N_hat. We also cannot go larger than rank_
            Ut = U(:,1:rank_);
            St = S(1:rank_,1:rank_);
            Vt = V(:,1:rank_);

            A_hat = St*Vt';
            A_dot = Ut'*A_dot';

            A_hat = A_hat';
            A_dot = A_dot';
        end
        %%

        % if options.rom.rom_type == "OpInf"
        if rom_type == "OpInf"
            [Diff,Conv] = OpInf_core(A_hat,A_dot);
        else
        % M = options.rom.M;
            M = size(A_dot,1);
            [operator_constraint,constraint_rhs] = get_constraints(M,rom_type);
            
            %% botch!!
            % operator_constraint = [];
            % constraint_rhs = [];
            %%

            A_dot_T = A_dot';

            % increments = [2];
            % O_2 = incremental_opinf(A_hat,A_dot_T,operator_constraint, constraint_rhs, increments);

            %% tryout: incremental OpInf
            % A_vec = sparse(kron(eye(M),A_hat'));
            % B_vec = A_dot_T(:);
            % 
            % NNhat = size(A_vec,2);
            % % inds1 = 1:2:NNhat;
            % % inds2 = 2:2:NNhat;
            % M1 = 1;
            % [inds1,inds2] = incremental_inds(M1,M);
            % 
            % O_s = zeros(NNhat,1);
            % 
            % 
            % % O_1 = lsqlin(A_vec(:,inds1),B_vec,[],[],operator_constraint(:,inds1),constraint_rhs);
            % % O_2 = lsqlin(A_vec(:,inds2),B_vec,[],[],operator_constraint(:,inds2),constraint_rhs);
            % % constraint-free
            % O_1 = lsqlin(A_vec(:,inds1),B_vec,[],[]);
            % O_s(inds1) = O_1;
            % O_2 = lsqlin(A_vec(:,inds2),B_vec - A_vec*O_s,[],[]);
            % 
            % O_c = lsqlin(A_vec,B_vec,[],[]);
            % 
            % O_s(inds2) = O_2;
            % 
            % norm(O_s-O_c)
            % 
            % norm(A_vec*O_c-B_vec)
            % norm(A_vec*O_s-B_vec)

            %% feasible basis approach (via linear system!)
            % A_vec = sparse(kron(eye(M),A_hat'));
            B_vec = A_dot_T(:);

            if rom_type == "EC-OpInf skew" && M>1
                T = sparse(feasible_basis(M));
                TT = blkdiag(eye(M^2), T);
                ordering = [reshape(1:M^2,M,M); M^2 + reshape(1:M^3,M^2,M)];
                TT = TT(ordering(:),:);
                % O_f = TT*((A_vec*TT)\B_vec);
                O_f = TT*((A_vecT(A_hat',TT))\B_vec);
                size(TT)

                O_ = O_f;
                
                % alternative to avoid rank deficiencies
                % note: this SVD could probably be combined with the
                % previous one, but SVD before multiplying TT is not
                % sufficient!
                % maybe multiple sequential SVDs are still useful as A_vec
                % caused a memory bottleneck before (and still does)
                % % [~,~,~,c_f2] = OpInf_SVD(full(A_vec*TT)',B_vec');
                % % O_f2 = TT*c_f2';

                % % norm(O_f - O_f2)

                %% debug stuff
                % O_o = lsqlin(sparse(kron(eye(M),(A_hat'))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
                % norm(O_f-O_o)
                % A_vec = sparse(kron(eye(M),A_hat'));
                % 
                % norm(A_vec*O_f - B_vec)
                % norm(A_vec*O_ - B_vec)
                %% 

                [Diff,Conv] = vec2ops(O_,M);
                

                cvx_begin
                    variable D(M,M) symmetric
                    minimize( norm(A_dot - [D,Conv]*A_hat) )
                    D <= 0
                cvx_end

                Diff = D;

            else
            



            %%

            % O_ = lsqlin(tall(sparse(kron(eye(M),(A_hat')))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            O_ = lsqlin(sparse(kron(eye(M),(A_hat'))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            % O_ = lsqlin(kron(eye(M),(A_hat')), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            % O_ = lsqlin(kron(eye(M),[A;A_kron]'), A_dot_T(:), [],[], [],[]); % lsqlin without constaints

                        [Diff,Conv] = vec2ops(O_,M);

            end


            %% botch: diffusion operator
            % A_hat_vec = kron(eye(M),(A_hat'));
            % A_dot_vec = A_dot_T(:);
            % disp("fmincon starting")
            % O_ = 0*O_; % zero initial solution
            % O_1 = fmincon(@(O) norm(A_hat_vec*symm_neg_def_diffusion(O,M) - A_dot_vec),O_, [],[], operator_constraint, constraint_rhs);
            % O_2 = symm_neg_def_diffusion(O_1,M);
            % [Diff,Conv] = vec2ops(O_2,M);
            %%
        end