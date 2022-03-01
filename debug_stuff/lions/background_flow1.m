function [V_delta,M_inbc,y_inbc] = background_flow1(options)

id_V_eff = find(options.grid.id_V_inbc);
id_p_eff = find(options.grid.id_p_inbc);

y_M = options.discretization.yM;
M_h = options.discretization.M;

y_inbc = -y_M(id_p_eff);
M_inbc = M_h(id_p_eff,id_V_eff);

sum(y_inbc)
% figure
% heatmap(M_inbc)

y_inbc(1) = 0;
% y_inbc(1) = 1;
% M_inbc(1,1) = 0;
% M_inbc(1,:) = 1; %nonzeros(options.grid.id_V_inbc_signed);
M_inbc(1,:) = nonzeros(options.grid.id_V_inbc_signed);


% figure
% heatmap(M_inbc)


V_delta_ = M_inbc\y_inbc;

NV = options.grid.NV;
V_delta = zeros(NV,1);
V_delta(id_V_eff) = V_delta_;

% vis_velo(V_delta,options)
% sum(V_delta)
sum(options.grid.id_V_inbc_signed.*V_delta)
17
%%
% y_inbc2 = y_inbc;
% y_inbc2(1) = 7/2;
% V_delta_2 = M_inbc\y_inbc2;
% 
% NV = options.grid.NV;
% V_delta2 = zeros(NV,1);
% V_delta2(id_V_eff) = V_delta_2;
% 
% vis_velo(V_delta2,options)
% sum(V_delta2)
% 
% %%
% y_inbc3 = y_inbc;
% M_inbc3 = M_inbc;
% 
% y_inbc3(1) = 0;
% % M_inbc3(1,:) = 1;
% M_inbc3(1,:) = nonzeros(options.grid.id_V_inbc_signed);
% 
% figure
% heatmap(M_inbc3)
% 
% 
% V_delta_3 = M_inbc3\y_inbc3;
% 
% NV = options.grid.NV;
% V_delta3 = zeros(NV,1);
% V_delta3(id_V_eff) = V_delta_3;
% 
% vis_velo(V_delta3,options)
% sum(V_delta3)
