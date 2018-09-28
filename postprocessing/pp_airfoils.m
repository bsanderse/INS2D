folder = 'results/airfoils/NACA2412_Re500/';
% file_list = {'ibm','surfaceCp','camberCp','chordCp','chordCL','chordCL_xfoilClCd'};
file_list = {'ibm','surfaceCp','chordCp'};

% color_list = {'k','r','b','m','c','g'};
color_list = {'k','r','r','r','b','b'};
dash_list = {'-','-','--','-.','-','--'};
linewidth = 1;
close all
for ii=1:length(file_list)
    
    clear uh vh p x_cl y_cl x_k y_k;
    load([folder file_list{ii} '.mat']);
    
    u = reshape(uh,Nux_in,Nuy_in);
    v = reshape(vh,Nvx_in,Nvy_in);
    pres = reshape(p,Npx,Npy);
    uwake1 = interp2(xin',yp,u',x_c+1.2,yp);
    uwake2 = interp2(xin',yp,u',x_c+3,yp);
    vwake1 = interp2(xp',yin,v',x_c+1.2,yin);
    vwake2 = interp2(xp',yin,v',x_c+3,yin);

    

    figure(2)
    plot(uwake1,yp,[color_list{ii} dash_list{ii}],'Linewidth',linewidth)
    hold on
    figure(3)
    plot(uwake2,yp,[color_list{ii} dash_list{ii}],'Linewidth',linewidth)
    hold on
    
    figure(4)
    plot(vwake1,yin,[color_list{ii} dash_list{ii}],'Linewidth',linewidth)
    hold on
    figure(5)
    hold on
    plot(vwake2,yin,[color_list{ii} dash_list{ii}],'Linewidth',linewidth)


    figure(6)
    list = 0:0.05:1.2;

    contour(xp,yp,qp',list)
    hold on
    if ( ii<=2 )
        patch(x_k,y_k,'w','Linewidth',2)      
%         plot(x_k,y_k,'k','Linewidth',2);
    else
        plot(x_cl,y_cl,'k','Linewidth',2);
    end
    axis equal
    axis([3.5 6.5 4 6])
    colorbar
    
    layout(6,14,'x','y');

    saveas(gcf,[folder file_list{ii} '_velocity.fig']);
    print(gcf,'-depsc',[folder file_list{ii} '_velocity.eps']);

    close(6)
    
    figure(7)
    list = -0.25:0.05:0.5;
    contour(xp,yp,pres',list);
    hold on
    if ( ii<=2 )
        patch(x_k,y_k,'w','Linewidth',2)
%         plot(x_k,y_k,'k','Linewidth',2);
    else
        plot(x_cl,y_cl,'k','Linewidth',2);
    end
    axis equal
    axis([3.5 6.5 4 6])
    colorbar
    
    layout(7,14,'x','y');

    saveas(gcf,[folder file_list{ii} '_pressure.fig']);
    print(gcf,'-depsc',[folder file_list{ii} '_pressure.eps']);    

    close(7)
    
end

%%
figure(2)
% legend('IBM','IBM','surface Cp','surface Cp','camber Cp','camber Cp','chord Cp','chord Cp','chord CL','chord CL')
% legend('IBM','surface Cp','camber Cp','chord Cp','chord CL','chord CL - Xfoil')
legend('method A','method B','method C','method D','method E','method F','Location','NorthWest')
layout(2,14,'u','y');
axis([0 1.1 4 6]);
grid on
box on
saveas(gcf,[folder 'uwake1.fig']);
print(gcf,'-depsc',[folder 'uwake1.eps']);  

figure(3)
% legend('IBM','surface Cp','camber Cp','chord Cp','chord CL','chord CL - Xfoil')
legend('method A','method B','method C','method D','method E','method F','Location','NorthWest')
layout(3,14,'u','y');
axis([0.5 1.1 4 6])
grid on
box on
saveas(gcf,[folder 'uwake2.fig']);
print(gcf,'-depsc',[folder 'uwake2.eps']);  

figure(4)
% legend('IBM','surface Cp','camber Cp','chord Cp','chord CL','chord CL - Xfoil')
legend('method A','method B','method C','method D','method E','method F','Location','NorthEast')
axis([-0.1 0.06 4 6])
layout(4,14,'v','y');
grid on
box on
saveas(gcf,[folder 'vwake1.fig']);
print(gcf,'-depsc',[folder 'vwake1.eps']);  

figure(5)
% legend('IBM','surface Cp','camber Cp','chord Cp','chord CL','chord CL - Xfoil')
legend('method A','method B','method C','method D','method E','method F','Location','SouthEast')
layout(5,14,'v','y');
grid on
box on
saveas(gcf,[folder 'vwake2.fig']);
print(gcf,'-depsc',[folder 'vwake2.eps']);  