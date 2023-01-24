function vis_p(p,options)

Np = options.grid.Np;

Npx = options.grid.Npx;
Npy = options.grid.Npy;

p_vis = reshape(p(1:Np),[Npx, Npy]);

figure
heatmap(p_vis')