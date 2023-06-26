                figure(117)
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx))/snapshots.k(1));
%                 plot(t_vec,k,'color',color,'displayname',name);

                if j==1 %Nsim
                    plot(t_vec,snapshots.k(snapshot_indx'),'color','k','linestyle','-','linewidth',linewidth,'displayname',"FOM");
                    hold on
                end

                plot(t_vec,k,'color',color,'linestyle',linestyle,'linewidth',linewidth,'displayname',name);
                hold on

%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1));
%                 semilogy(t_vec,abs(k-k(1))/k(1));

%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1),color,'displayname',"ROM M="+M+suffix);
%                 semilogy(t_vec,abs(k-k(1))/k(1),color,'displayname',"ROM M="+M+suffix);

%                 hold on
%                 set(gca,'Yscale','log')
%                 ylabel('energy error (K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)');
%                 ylabel('energy error (K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)');
                ylabel('kinetic energy')

%                 legend('(K_{ROM}(t)-K_{FOM}(t))/K_{FOM}(0)','(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)','(K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)')
%                 legend('(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)')
%                 title('error in kinetic energy ROM');
                  legend('show')

                      set(gcf, 'Position', [100, 100, 400, 300])
