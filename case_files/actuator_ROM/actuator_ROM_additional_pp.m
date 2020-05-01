% additional postprocessing, run after post_processing.m is completed in
% order to get certain variables

%% script to plot Vbc function

% Vbc = options.rom.Vbc;

% u-component
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);
[~,c]=contour(options.grid.xin,options.grid.yp,reshape(Vbc(1:options.grid.Nu),options.grid.Nux_in,options.grid.Nuy_in)');
c.LineWidth=2;
axis equal
axis([x1 x2 y1 y2]);
colorbar
xlabel('x')
ylabel('y')

% export_fig('actuator_uBC_extended','-pdf','-transparent')

% v-component
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);
[~,c]=contour(options.grid.xp,options.grid.yin,reshape(Vbc(options.grid.Nu+1:end),options.grid.Nvx_in,options.grid.Nvy_in)',20);
c.LineWidth=2;
axis equal
axis([x1 x2 y1 y2]);
colorbar
xlabel('x')
ylabel('y')

% export_fig('actuator_vBC_extended','-pdf','-transparent')

%% script to plot final pressure and velocity
[up,vp,qp] = get_velocity(V,t,options);
list = linspace(0.5,1.15,20);
% list = 20;
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);
% pcolor(xp,yp,qp')
[~,c]=contour(xp,yp,qp',list);
c.LineWidth = 2;
axis equal
axis([x1 x2 y1 y2]);
colorbar
xlabel('x')
ylabel('y')

% export_fig('ROM_velocity_field_M10_t20','-pdf','-transparent')
% export_fig('FOM_velocity_field_t20','-pdf','-transparent')

Npx = options.grid.Npx;
Npy = options.grid.Npy;    
pres = reshape(p,Npx,Npy);
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);
list = linspace(-0.1,0.1,20);
[~,c]=contour(xp,yp,pres',list);
c.LineWidth = 2;
axis equal
axis([x1 x2 y1 y2]);
colorbar
xlabel('x')
ylabel('y')
% export_fig('ROM_pressure_field_M10_t20','-pdf','-transparent')
% export_fig('FOM_pressure_field_t20','-pdf','-transparent')

%% script to plot mass conservation 
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);

plot(t_vec,snapshots.maxdiv(snapshot_indx),'LineWidth',2);
hold on
plot(t_vec,maxdiv,'LineWidth',2);
xlabel('t')
ylabel('maximum divergence of velocity field');
grid
legend('FOM','ROM');
set(gca, 'YScale', 'log')
set(gca, 'FontSize',14)

export_fig('error_div_actuator_unsteady_M10','-pdf','-transparent')

%% script to plot errors
figure
set(gcf,'color','w');
set(gca,'LineWidth',1);
plot(t_vec,error_V_2,'LineWidth',2);
hold on
plot(t_vec,error_V_best_2,'LineWidth',2);
legend('Velocity - ROM','Velocity - best approximation');

if (options.rom.pressure_recovery == 1)
    plot(t_vec,error_p_2,'LineWidth',2);
    hold on
    plot(t_vec,error_p_best_2,'LineWidth',2);
%     legend('L_{inf} error in ROM velocity','L_2 error in ROM velocity','L_{inf} error in ROM pressure')
    legend('Velocity - ROM','Velocity - best approximation','Pressure - ROM','Pressure - best approximation')
%     legend('Velocity','Pressure')
end
set(gca, 'YScale', 'log')
set(gca, 'FontSize',14)
grid
xlabel('t')
ylabel('L2 error')
export_fig('error_velocity_actuator_unsteady','-pdf','-transparent')
