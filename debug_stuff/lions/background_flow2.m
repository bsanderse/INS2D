function [V_delta,matrix,y_inbc] = background_flow2(options)

id_V_eff = find(options.grid.id_V_inbc2);
id_p_eff = find(options.grid.id_p_inbc2);

y_M = options.discretization.yM;
M_h = options.discretization.M;

y_inbc = -y_M(id_p_eff); sum(y_inbc)
M_inbc = M_h(id_p_eff,id_V_eff);

G_inbc = -M_inbc';

Om_inv = options.grid.Om_inv;

Om_inv_inbc = Om_inv(id_V_eff);

L_inbc = M_inbc*(Om_inv_inbc.*G_inbc);

%% make problem well-posed  -> apparently not required!
% L_inbc(1,:) = 1;
% y_inbc(1) = 0;
%%

sum(y_inbc)

L_inbc_inv = inv(L_inbc);
matrix = Om_inv_inbc.*(G_inbc*L_inbc_inv);
V_delta_ = Om_inv_inbc.*(G_inbc*(L_inbc\y_inbc));

%% debugging
% norm(M_inbc*V_delta_-y_inbc)
% A =M_inbc*matrix;
% figure
% heatmap(A)


%%

NV = options.grid.NV;

V_delta = zeros(NV,1);
V_delta(id_V_eff) = V_delta_;

% % figure
% % heatmap(M_inbc)
% 
% y_inbc(1) = 0;
% % y_inbc(1) = 1;
% % M_inbc(1,1) = 0;
% % M_inbc(1,:) = 1; %nonzeros(options.grid.id_V_inbc_signed);
% M_inbc(1,:) = nonzeros(options.grid.id_V_inbc_signed);
% 
% % figure
% % heatmap(M_inbc)
% 
% %%
% 
% % V_delta_ = M_inbc\y_inbc;
% % 
% % NV = options.grid.NV;
% % V_delta = zeros(NV,1);
% % V_delta(id_V_eff) = V_delta_;
% % 
% % % vis_velo(V_delta,options)
% % % sum(V_delta)
% % sum(options.grid.id_V_inbc_signed.*V_delta)
% % 17