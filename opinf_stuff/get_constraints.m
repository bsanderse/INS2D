function [operator_constraint,constraint_rhs] = get_constraints(M,rom_type)
            convection_constraint = [];
            % switch options.rom.rom_type
            switch rom_type
                case "EC-OpInf Koike"
                    convection_constraint = [convection_constraint; three_term_prop_constraint(M)];
                    convection_constraint = [convection_constraint; ambiguity_constraint(M)];
                case "EC-OpInf skew"
                    convection_constraint = [convection_constraint; block_skewsymm_constraint(M)'];
                    convection_constraint = [convection_constraint; triple_ambiguity_constraint(M)'];
            end
            n_constr = size(convection_constraint,1);
            % operator_constraint_ = [zeros(n_constr,M^2) convection_constraint];
            operator_constraint_ = [sparse(n_constr,M^2) convection_constraint];
            ordering = [reshape(1:M^2,M,M); M^2 + reshape(1:M^3,M^2,M)];
            operator_constraint = operator_constraint_(:,ordering(:));

            %% remove all zero rows
            operator_constraint = operator_constraint(sum(abs(operator_constraint),2)~=0,:);
            n_constr = size(operator_constraint,1);
            %%

            constraint_rhs = sparse(n_constr,1);