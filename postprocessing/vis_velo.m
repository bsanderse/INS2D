function vis_velo(V,options)

Nu = options.grid.Nu;
Nv = options.grid.Nv;

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;

Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

u_vis = reshape(V(1:Nu),[Nux_in, Nuy_in]);
v_vis = reshape(V(Nu+1:end),[Nvx_in, Nvy_in]);

figure
subplot(2,1,1)
heatmap(u_vis')

subplot(2,1,2)
heatmap(v_vis')