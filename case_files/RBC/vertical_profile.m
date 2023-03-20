Tmean_temp = reshape(Tmean,Npx,Npy);
sum=zeros(Ny,1);
for i=1:Ny
    for j=1:Nx
        sum(i)=sum(i)+Tmean_temp(j,i);
    end
    sum(i)=sum(i)/Nx;
end

figure(201)
plot(sum,yp,'r');
grid on
hold on
xlim([0 1]);
ylim([0 1]);
xlabel('mean temperature');
ylabel('y');
set(gcf,'color','w');
set(gca,'LineWidth',1,'FontSize',14);
