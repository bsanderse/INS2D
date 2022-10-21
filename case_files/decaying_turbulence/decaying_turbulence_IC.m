function [u,v,p,options] = decaying_turbulence_IC(t,options)
% initial velocity field

yin = options.grid.yin;
xin = options.grid.xin;
Npx = options.grid.Npx;
Npy = options.grid.Npy;

psi = zeros(Npx,Npy);
nv = 16;

dx = 0.01*randn(nv,nv);
dy = 0.01*randn(nv,nv);
if (options.output.save_results==1)
    writematrix(dx,[options.output.path_results '/dx.txt'])
    writematrix(dy,[options.output.path_results '/dy.txt'])
end
xpsi = kron(ones(Npy,1),xin);
ypsi = kron(yin,ones(Npx,1));
xpsi = reshape(xpsi,Npx,Npy);
ypsi = reshape(ypsi,Npx,Npy);

for j = 1:nv
    for k =1:nv
        
        psi = psi + 5e-2 * ((-1)^(k+j)) * ...
            exp(-2000*( (xpsi + dx(k,j) - k/(nv+1)).^2 + (ypsi + dy(k,j) - j/(nv+1)).^2 ) );
        
    end
end


V   = get_velocity_from_streamfunction(psi(:),t,options);
u   = V(options.grid.indu);
v   = V(options.grid.indv);
p   = zeros(Npx,Npy);

    
end