                figure(113)
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx))/snapshots.k(1));
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx')),'color',color,'displayname',name); ylabel('energy error')
                k_avg = sum(k)/numel(k)
                semilogy(t_vec,abs(k - snapshots.k(snapshot_indx'))/k_avg,'color',color,'displayname',name); ylabel('relative energy error')

%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1));
%                 semilogy(t_vec,abs(k-k(1))/k(1));

%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1),color,'displayname',"ROM M="+M+suffix);
%                 semilogy(t_vec,abs(k-k(1))/k(1),color,'displayname',"ROM M="+M+suffix);

                hold on
%                 set(gca,'Yscale','log')
%                 ylabel('energy error (K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)');
%                 ylabel('energy error (K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)');
%                 ylabel('energy error')

%                 legend('(K_{ROM}(t)-K_{FOM}(t))/K_{FOM}(0)','(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)','(K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)')
%                 legend('(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)')
%                 title('error in kinetic energy ROM');
                  legend('show')

                      set(gcf, 'Position', [100, 100, 400, 300])
