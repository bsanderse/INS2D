M_list = [2 2 4 4 8 8 16 16 32 32];

svd = [   1.838781550000000   1.845571699000000   1.865703083000000   1.783140290000000  1.740951802000000   1.815695423000001  1.900404364000000   1.855649407000000  1.749046473000000   1.983229209000000];

precompute =  [0.251733964000000   0.248007699000000    0.465944821000000   0.464967988000000   1.286483787000000   1.248063014000000   4.100600315000000   4.040544959999999   16.991264284000003  15.634421606000000];

time_loop =   [ 2.364090564000001   1.830272943000001    1.880433734000000   1.860130696000000   1.924760972000000   1.900528199999999    2.122761275000000   2.145150909000000   2.754535253000000   2.612073565000003];
% FOM has been run 3 times:
% results are for the time_loop
fom = [     1.004640919380000e+02 1.082086400460000e+02   1.106526816070000e+02];


avg =@(list) 0.5*(list(1:2:end) + list(2:2:end));

figure
set(gca,'LineWidth',1);
M_plot = M_list(1:2:end);
svd_plot = avg(svd);

precompute_plot = avg(precompute);

time_plot = avg(time_loop);

loglog(M_plot,svd_plot,'s-','LineWidth',2);
hold on
loglog(M_plot,precompute_plot,'s-','LineWidth',2);
loglog(M_plot,time_plot,'s-','LineWidth',2);
loglog([1;60],[mean(fom);mean(fom)],'--','LineWidth',2);

grid
legend('Offline - SVD','Offline - Precomputing operators','Online','FOM','Location','southeast');
xlabel('M');
ylabel('CPU time [s]');
xlim([1 60]);
ylim([0.5e-1 3e2]);
set(gca,'FontSize',14);

%%
export_fig('CPUtimings_shearlayer_unsteady','-pdf','-transparent')
