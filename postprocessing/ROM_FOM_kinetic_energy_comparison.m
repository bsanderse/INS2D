% ROM_FOM_kinetic_energy_comparison

ts = t_start:dt:t_end;
k_FOM = snapshots.k;

figure(53)
if j==1
    plot(ts,k_FOM,'displayname','k FOM')
    hold on
end

plot(ts,k,'displayname',"k ROM M="+M+suffix)
legend('show')

%%
figure(54)
% plot(ts,k_FOM-k,'displayname',"k FOM - k ROM M="+M+suffix)
semilogy(ts,abs(k_FOM-k),'displayname',"k FOM - k ROM M="+M+suffix)
if j==1
    hold on
end

legend('show')

% 17

