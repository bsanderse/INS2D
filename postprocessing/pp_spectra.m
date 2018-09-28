folder_results = 'results/regularization/nomodel_Re500/';

fac = Npx/10;
detailed_plot = 0;

%%
figure
surf(10*dt:10*dt:t_end,freq,Ek1)
set(gca,'ZScale','log','YScale','log')
logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ek1_N',num2str(Nx),'_t',num2str(t_end)],'fig');

%%
if (detailed_plot==1)
figure
surf(10*dt:10*dt:t_end,freq(1:end/fac),Ek1(1:end/fac,:))
% set(gca,'ZScale','log','YScale','log')
% logzplot colorbar
shading interp
caxis([1e-2 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ek1_N',num2str(Nx),'_t',num2str(t_end),'_detail'],'fig');
end

%%
figure
surf(10*dt:10*dt:t_end,freq,Ek2)
set(gca,'ZScale','log','YScale','log')
logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ek2_N',num2str(Nx),'_t',num2str(t_end)],'fig');

%%
if (detailed_plot==1)

figure
surf(10*dt:10*dt:t_end,freq(1:end/fac),Ek2(1:end/fac,:))
% set(gca,'ZScale','log','YScale','log')
% logzplot colorbar
shading interp
caxis([1e-2 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ek2_N',num2str(Nx),'_t',num2str(t_end),'_detail'],'fig');

end

%%
figure
surf(10*dt:10*dt:t_end,freq,Ew1)
set(gca,'ZScale','log','YScale','log')
logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ew1_N',num2str(Nx),'_t',num2str(t_end)],'fig');
    

%%
if (detailed_plot==1)

figure
surf(10*dt:10*dt:t_end,freq(1:end/fac),Ew1(1:end/fac,:))
% set(gca,'ZScale','log','YScale','log')
% logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ew1_N',num2str(Nx),'_t',num2str(t_end),'_detail'],'fig');
    
end

%%
figure
surf(10*dt:10*dt:t_end,freq,Ew2)
set(gca,'ZScale','log','YScale','log')
logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ew2_N',num2str(Nx),'_t',num2str(t_end)],'fig');

%%
if (detailed_plot==1)

figure
surf(10*dt:10*dt:t_end,freq(1:end/fac),Ew2(1:end/fac,:))
% set(gca,'ZScale','log','YScale','log')
% logzplot colorbar
shading interp
caxis([1e-5 1e0])
% view([140 30]);
view(2)
saveas(gcf,[folder_results,'Ew2_N',num2str(Nx),'_t',num2str(t_end),'_detail'],'fig');
end