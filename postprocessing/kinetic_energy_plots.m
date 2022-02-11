% kinetic_energy_plots
figure
ts = t_start:dt:t_end;
delta_k = (k(2:end)-k(1:end-1))/dt;
plot(ts(1:end-1),delta_k,'displayname','k(n+1)-k(n)/dt');
hold on
plot(ts(1:end-1),k_diff(1:end-1),'displayname','k diff');
plot(ts(1:end-1),k_conv(1:end-1),'displayname','k conv');
plot(ts(1:end-1),k_pres(1:end-1),'displayname','k pres');
plot(ts(1:end-1),k_force(1:end-1),'displayname','k force');
plot(ts(1:end-1),k_diffBC(1:end-1),'displayname','k diff BC');
plot(ts(1:end-1),k_presBC(1:end-1),'displayname','k pres BC');

k_sum = k_diff + k_conv + k_pres + k_force + k_diffBC + k_presBC;
plot(ts(1:end-1),k_sum(1:end-1),'displayname','k sum');


legend('show')

figure
plot(ts(1:end-1),delta_k'-k_sum(1:end-1),'displayname','k(n+1)-k(n)/dt - k sum');
hold on
plot(ts(1:end-1),delta_k-k_sum2(1:end-1),'displayname','k(n+1)-k(n)/dt - k sum');
legend('show')

figure
semilogy(ts(1:end-1),abs(delta_k'-k_sum(1:end-1)),'displayname','k(n+1)-k(n)/dt - k sum');
hold on
semilogy(ts(1:end-1),abs(delta_k-k_sum2(1:end-1)),'displayname','k(n+1)-k(n)/dt - k sum');
legend('show')
17