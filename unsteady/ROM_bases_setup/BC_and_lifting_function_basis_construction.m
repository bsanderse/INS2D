%% compute boundary condition approximation and inhomogeneous ROM basis
if (options.rom.bc_recon == 3) || (options.rom.bc_recon == 5) 
        if options.rom.rom_bc == 2
            [U_bc,S_bc,V_bc] = svd(X_bc,'econ');
            if (options.visualization.show_sigmas == 1)
                Sigma_bc = diag(S_bc);
                figure
                semilogy(Sigma_bc/Sigma_bc(1),'s','displayname', 'singular values Vbc snapshot matrix');
                ylabel("\sigma_i/\sigma_1")
                xlabel("mode index")
                title('singular values')
                legend('show')
                if (exist('fig_destination') && j==Nsim)
                    savefig([fig_destination '/singular values'])
                end
            end
        else
            U_bc = X_bc/norm(X_bc);
        end
        if options.verbosity.equivalence_cheat ~= 2
            cond_fac = 10^-6;
            X_bc_rank = sum(abs(Sigma_bc/Sigma_bc(1))>cond_fac);
            if Mbc > X_bc_rank
                Mbc = X_bc_rank;
            end

            phi_bc = U_bc(:,1:Mbc);
            options.rom.phi_bc = phi_bc;
        end

        if (options.rom.bc_recon == 3) || options.verbosity.equivalence_cheat == 1
            for jj = 1:Mbc
                yBC = phi_bc(:,jj);
                Y_M(:,jj) = get_yM(options,yBC);
            end
            L = options.discretization.A;
            Gx   = options.discretization.Gx;
            Gy   = options.discretization.Gy;
            G = [Gx;Gy];
            Om = options.grid.Om;
            Om_inv = options.grid.Om_inv;
            %     tilde_phi_inhom = Om_inv.*(G*(L\Y_M));
            tilde_phi_inhom = -Om_inv.*(G*(L\Y_M)); %pfusch

            [Q_inhom,R_inhom] = qr(sqrt(Om).*tilde_phi_inhom,0); %alternative: take first vec of tilde phi inhom
            M_inhom = rank(tilde_phi_inhom);
            Q_1_inhom = -Q_inhom(:,1:M_inhom);
            R_inhom = -R_inhom(1:M_inhom,:);
            phi_inhom = sqrt(Om_inv).*Q_1_inhom;

            options.rom.phi_inhom = phi_inhom;
            options.rom.R_inhom = R_inhom;
%%
%             warning('phi_inhom computation manipulated')
%             M_h = options.discretization.M;
%             [Qq,Rr] = qr((M_h*phi_inhom)');
%             rankk = rank(M_h*phi_inhom);
%             Qq1 = Qq(:,1:rankk);
%             phi_inhom_star = phi_inhom*Qq1;
%             R_inhom_star = Qq1'*R_inhom;
% 
%             options.rom.phi_inhom = phi_inhom_star;
%             options.rom.R_inhom = R_inhom_star;
%%
        else
            options.rom.phi_inhom = 0;
            options.rom.R_inhom = 0;
        end
end

if options.verbosity.equivalence_cheat == 1
    B = [B phi_inhom];
    options.rom.B = B;
    M = M + M_inhom;
    options.rom.M = M;

    %%
    if options.rom.bc_recon == 5
        % for bc_recon == 5, existence of a solution to the ROM PPE requires in
        % general the invertibility of \hat L = \hat M \hat G <=>  \hat M has
        % full rank
        if M < Mp 
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = M;
        end
        M_h = options.discretization.M;
        Bp = Wp(:,1:Mp);
%         while rank(Bp'*M_h*B)<Mp % prone to machine precision problems
%         rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
%         cond_fac = 10^-8;
        cond_fac = 10^-6;
        sing_vals = svd(Bp'*M_h*B);
        rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        while rank_<Mp
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = rank_;
            Bp = Wp(:,1:Mp);
%             rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
            sing_vals = svd(Bp'*M_h*B);
            rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        end
    else
        Bp = Wp(:,1:Mp);
    end 
    options.rom.Bp = Bp;
    %%
end