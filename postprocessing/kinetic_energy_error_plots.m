                figure(113)
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx))/snapshots.k(1));
                semilogy(t_vec,abs(k - snapshots.k(snapshot_indx')),color,'displayname',"ROM M="+M+suffix);
%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1));
%                 semilogy(t_vec,abs(k-k(1))/k(1));

%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1),color,'displayname',"ROM M="+M+suffix);
%                 semilogy(t_vec,abs(k-k(1))/k(1),color,'displayname',"ROM M="+M+suffix);

                hold on
%                 set(gca,'Yscale','log')
%                 ylabel('energy error (K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)');
%                 ylabel('energy error (K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)');
                ylabel('energy error')

%                 legend('(K_{ROM}(t)-K_{FOM}(t))/K_{FOM}(0)','(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)','(K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)')
%                 legend('(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)')
%                 title('error in kinetic energy ROM');
                  legend('show')