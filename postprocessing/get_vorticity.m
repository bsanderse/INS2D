function [omega] = get_vorticity(V,t,options)
% vorticity values at pressure midpoints
% this should be consistent with operator_postprocessing

BC  = options.BC;

Nu  = options.grid.Nu;
Nv  = options.grid.Nv;
Nux_in  = options.grid.Nux_in;
Nvy_in  = options.grid.Nvy_in;

Nx  = options.grid.Nx;
Ny  = options.grid.Ny;

Wv_vx = options.discretization.Wv_vx;
Wu_uy = options.discretization.Wu_uy;

uh = V(1:Nu);
vh = V(Nu+1:Nu+Nv);


if (strcmp(BC.u.left,'per') && strcmp(BC.v.low,'per'))
    uh_in = uh;
    vh_in = vh;
else
    % velocity at inner points
    diagpos = 0;
    if (strcmp(BC.u.left,'pres') )
        diagpos = 1;
    end
    if (strcmp(BC.u.right,'per') && strcmp(BC.u.left,'per') ) % like pressure left
        diagpos = 1;
    end
    
    B1D    = spdiags(ones(Nx-1,1),diagpos,Nx-1,Nux_in);
    B2D    = kron(speye(Ny),B1D);
    
    uh_in  = B2D*uh;
    
    
    diagpos = 0;
    if (strcmp(BC.v.low,'pres') )
        diagpos = 1;
    end
    if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per') ) % like pressure low
        diagpos = 1;
    end
    
    B1D    = spdiags(ones(Ny-1,1),diagpos,Ny-1,Nvy_in);
    B2D    = kron(B1D,speye(Nx));
    
    vh_in  = B2D*vh;
end

omega  = Wv_vx*vh_in - Wu_uy*uh_in;

end