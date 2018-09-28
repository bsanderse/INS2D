% plot(t_plot,k(1:n),'kx-');
% hold on
% % plot(t_plot(1:n),k_error,'kx-');
% plot(t_plot(1:n),k_exact(1:n),'ro-')

% 
% xlabel('$t$','Interpreter','Latex');
% ylabel('$\epsilon_k$','Interpreter','Latex');
% legend('kinetic energy');
% % set(gca,'XTick',0:2:16);
% % set(gca,'XTickLabel',{'0','2','4','6','8','6','4','2','0'})
% grid
% k_error_list(j) = abs(k(end)-k(1))


% % figure
% % t_plot = [0:abs(dt):t_end/2 t_end/2:dt:0];
t_plot = 0:abs(dt):t_end;
% % t_plot(n) = t; 
% % t_plot = 1:n;
% % plot(t_plot,umom(1:n)-umom(1),'bx-');
% % hold on
% % plot(t_plot,vmom(1:n)-vmom(1),'rx-');
% % legend('u-momentum','v-momentum');
% % 
% figure
k_error = (k(1:n)-k(1))/k(1);


% keyboard;
% % time reversibility error
% % uh_end = uh;
% % vh_end = vh;
% % initialize;
% % uh_start = u(:);
% % vh_start = v(:);
% % 
% % rev_error(end) = sum(Omu.*(uh_end - uh_start).^2) + sum(Omv.*(vh_end - vh_start).^2)
% 
%

colors1 = {'kx-','rx-','bx-','gx-','mx-','yx-','cx-'};
colors2 = {'ko--','ro--','bo--','go--','mo--','yo--','co--'};


figure(2)
% plot(t_plot(nt/2+2:end),rev_error(nt/2+2:end),'kx-');
plot(t_plot(1:n),k_error(1:n),char(colors1(jj)))
hold on

% xlabel('$t$','Interpreter','Latex');
% ylabel('$\epsilon_t$','Interpreter','Latex');
% legend('time reversibility');
set(gca,'XTick',0:2:16);
set(gca,'XTickLabel',{'0','2','4','6','8','6','4','2','0'})
% grid
% rev_error_list(j) = abs(rev_error(end))
% 
% %%
figure(3)
enstrophy_error = (enstrophy(1:n)-enstrophy(1))/enstrophy(1);
vorticity_error = omega_total;
plot(t_plot(1:n),enstrophy_error,char(colors1(jj)));
hold on
% plot(t_plot(1:n),omega_total,char(colors2(j)));
grid
enstrophy_error_end(j,jj) = abs(enstrophy_error(end))
enstrophy_error_mid(j,jj) = abs(enstrophy_error(nt/2+1))
% legend('enstrophy error','vorticity error')

if (exist('nonlinear_its') && length(nonlinear_its(2:end))==length(t_plot(2:n)))
figure(4)
plot(t_plot(2:n),nonlinear_its(2:end),char(colors1(jj)));
hold on
end

reversibility_error(j,jj) = max(abs([uh;vh]-[uh_start(:);vh_start(:)]))
% energy_error(j,jj) = k_error(nt/2+1)
energy_error(j,jj) = k_error(end)