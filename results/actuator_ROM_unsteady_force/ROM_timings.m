M_list = [5 5 10 10 20 20 40 40 80 80];

svd = [   0.363540362000000; 0.400579872000000; 0.369285986000000; 0.388155839000000; 0.377476190000000;  0.368451502000000;  0.385662822000000; 0.396618829000000; 0.371343766000000;  0.382594722000000];

precompute = [0.407085773000000; 0.392370833000000; 1.075285911000000; 1.089639901000000; 3.411268463000000; 3.471748365999999;  14.729502216000000; 14.014026250000001; 63.097201753999997; 60.774031010000002];

time_loop = [0.490072280000000;  0.511090747000000; 0.558533572000000; 0.558945277000000; 0.687801013000000;  0.769460111000000; 1.033060843999998; 1.225009048000000;  2.170482172999996; 2.386989877000005];

% FOM has been run 3 times:
fom = [73.664876351999993; 76.077162514999998; 81.306319036999994];


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
loglog([2;100],[mean(fom);mean(fom)],'--','LineWidth',2);

grid
legend('Offline - SVD','Offline - Precomputing operators','Online','FOM','Location','west');
xlabel('M');
ylabel('CPU time [s]');
xlim([2 100]);
set(gca,'FontSize',14);

%%
export_fig('CPUtimings_actuator_unsteady','-pdf','-transparent')
