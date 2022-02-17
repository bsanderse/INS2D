% kinetic_energy_plots

ts = t_start:dt:t_end;
% delta_k = (k(2:end)-k(1:end-1))/dt;
delta_k = (k(2:end)-k(1:end-1));

k_sum = k_diff + k_conv + k_pres + k_force + k_diffBC + k_presBC + k_obc;

% k_sum = k_analysis.k_diff + k_analysis.k_conv ...
%       + k_analysis.k_pres + k_analysis.k_force ...
%       + k_analysis.k_diffBC + k_analysis.k_presBC;

%% components
figure
plot(ts(1:end-1),delta_k,'x','displayname','k(n+1)-k(n)/dt');
hold on
plot(ts(1:end-1),k_diff(2:end),'displayname','k diff');
plot(ts(1:end-1),k_conv(2:end),'displayname','k conv');
plot(ts(1:end-1),k_pres(2:end),'displayname','k pres');
plot(ts(1:end-1),k_force(2:end),'displayname','k force');
plot(ts(1:end-1),k_diffBC(2:end),'displayname','k diff BC');
plot(ts(1:end-1),k_presBC(2:end),'displayname','k pres BC');
plot(ts(1:end-1),k_obc(2:end),'displayname','k pres BC');

plot(ts(1:end-1),k_sum(2:end),'s','displayname','k sum');

if options.rom.rom == 1
    snapshots.delta_k = (snapshots.k(2:end)-snapshots.k(1:end-1));
    
    plot(ts(1:end-1),snapshots.delta_k,'x','displayname','FOM k(n+1)-k(n)/dt');
    hold on
    plot(ts(1:end-1),snapshots.k_diff(2:end),'displayname','FOM k diff');
    plot(ts(1:end-1),snapshots.k_conv(2:end),'displayname','FOM k conv');
    plot(ts(1:end-1),snapshots.k_pres(2:end),'displayname','FOM k pres');
    plot(ts(1:end-1),snapshots.k_force(2:end),'displayname','FOM k force');
    plot(ts(1:end-1),snapshots.k_diffBC(2:end),'displayname','FOM k diff BC');
    plot(ts(1:end-1),snapshots.k_presBC(2:end),'displayname','FOM k pres BC');
    plot(ts(1:end-1),snapshots.k_obc(2:end),'displayname','FOM k pres BC');
 
    snapshots.k_sum = snapshots.k_diff + snapshots.k_conv ...
        + snapshots.k_pres + snapshots.k_force + snapshots.k_diffBC ...
        + snapshots.k_presBC + snapshots.k_obc;
    
    plot(ts(1:end-1),snapshots.k_sum(2:end),'s','displayname','k sum');
end

legend('show')
%%
% 
% figure
% plot(ts(1:end-1),delta_k'-k_sum(1:end-1),'displayname','k(n+1)-k(n)/dt - k sum');
% hold on
% plot(ts(1:end-1),delta_k-k_sum2(1:end-1),'displayname','k(n+1)-k(n)/dt - k sum');
% legend('show')

figure

semilogy(ts(1:end-1),abs(delta_k'-k_sum(2:end)),'displayname','k(n+1)-k(n)/dt - k sum');
hold on
semilogy(ts(1:end-1),abs(delta_k-k_sum2/2),'displayname','k(n+1)-k(n)/dt - k sum2');
legend('show')
17

% abs(delta_k'-k_sum(2:end)*dt) + abs(delta_k-k_sum2/2)
% semilogy(ts(1:end-1),abs(delta_k'-k_sum(2:end)*dt)' + abs(delta_k-k_sum2/2),'displayname','fun')

%% components
% figure
% plot(ts(1:end-1),delta_k,'displayname','k(n+1)-k(n)/dt');
% hold on
% plot(ts(1:end-1),k_analysis.k_diff(1:end-1),'displayname','k diff');
% plot(ts(1:end-1),k_analysis.k_conv(1:end-1),'displayname','k conv');
% plot(ts(1:end-1),k_analysis.k_pres(1:end-1),'displayname','k pres');
% plot(ts(1:end-1),k_analysis.k_force(1:end-1),'displayname','k force');
% plot(ts(1:end-1),k_analysis.k_diffBC(1:end-1),'displayname','k diff BC');
% plot(ts(1:end-1),k_analysis.k_presBC(1:end-1),'displayname','k pres BC');
% 
% plot(ts(1:end-1),k_sum(1:end-1),'displayname','k sum');


17

