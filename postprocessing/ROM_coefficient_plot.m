

figure(567)
for k = 1:options.rom.M
plot(t_vec,R_total(:,k),'displayname',num2str(k));
hold on
end

xlabel('t')
ylabel('ROM coefficients')
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show')


figure(568)
for k = 1:options.rom.Mbc
plot(t_vec,abc_total(:,k),'displayname',num2str(k));
hold on
end

xlabel('t')
ylabel('BC approximation coefficients')
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show')

figure(569)
for k = 20:21
plot(t_vec,abc_total(:,k),'displayname',num2str(k));
hold on
end

xlabel('t')
ylabel('BC approximation coefficients')
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show')

figure(47887)
% semilogy(sum(abs(abc_total),1))
% semilogy(sum(abs(abc_total),1)/numel(t_vec))
% sums = sum(abs(abc_total),1);

sums = sum(abs(abc_total),1);
% sums = sum(abc_total.^2,1);
semilogy(sums/sums(1),'displayname', "a_{bc}/a_{bc}(1) time average")

sums = sum(abc_total.^2,1);
semilogy(sums/sums(1),'displayname', "(a_{bc}/a_{bc}(1))^2 time average")

