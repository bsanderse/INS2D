function [Diff, Conv] = OpInf(A_dot,A_hat,rom_type) 

    %% botch
        if true
            [U,S,V] = svd(A_hat',"econ");

            rank = sum(abs(diag(S))>sqrt(eps))
            Ut = U(:,1:rank);
            St = S(1:rank,1:rank);
            Vt = V(:,1:rank);

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
            convection_constraint = [];
            % switch options.rom.rom_type
            switch rom_type
                case "EC-OpInf Koike"
                    convection_constraint = [convection_constraint; three_term_prop_constraint(M)];
                    convection_constraint = [convection_constraint; ambiguity_constraint(M)];
                case "EC-OpInf skew"
                    convection_constraint = [convection_constraint; block_skewsymm_constraint(M)];
                    convection_constraint = [convection_constraint; triple_ambiguity_constraint(M)];
            end
            n_constr = size(convection_constraint,1);
            operator_constraint_ = [zeros(n_constr,M^2) convection_constraint];
            ordering = [reshape(1:M^2,M,M); M^2 + reshape(1:M^3,M^2,M)];
            operator_constraint = operator_constraint_(:,ordering(:));
            constraint_rhs = zeros(n_constr,1);

            A_dot_T = A_dot';
            % O_ = lsqlin(tall(sparse(kron(eye(M),(A_hat')))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            O_ = lsqlin(sparse(kron(eye(M),(A_hat'))), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            % O_ = lsqlin(kron(eye(M),(A_hat')), A_dot_T(:), [],[], operator_constraint, constraint_rhs);
            % O_ = lsqlin(kron(eye(M),[A;A_kron]'), A_dot_T(:), [],[], [],[]); % lsqlin without constaints
            [Diff,Conv] = vec2ops(O_,M);

            %% botch: diffusion operator
            % A_hat_vec = kron(eye(M),(A_hat'));
            % A_dot_vec = A_dot_T(:);
            % O_1 = fmincon(@(O) norm(A_hat_vec*symm_neg_def_diffusion(O,M) - A_dot_vec),O_, [],[], operator_constraint, constraint_rhs);
            % O_2 = symm_neg_def_diffusion(O_1,M);
            % [Diff,Conv] = vec2ops(O_2,M);
            %%
        end