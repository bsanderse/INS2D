figure(5)

t_vec = t_start:dt:t_end;


K_h0 = snapshots.k(1); % K_h(0) initial FOM energy

semilogy(t_vec,abs(k - K_h0)/K_h0)
hold on

legend('show','NumColumns',3,'Orientation','vertical')


set(gcf, 'Position', [100, 100, 400, 300])
grid on
ylabel('relative energy error')
xlabel('t')
