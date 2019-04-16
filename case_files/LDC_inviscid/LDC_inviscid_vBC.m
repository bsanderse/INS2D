function v = LDC_inviscid_vBC(x,y,t,options)
% boundary conditions for u for LDC

    x1 = options.grid.x1;
    x2 = options.grid.x2;

    v = zeros(length(x)*length(y),1);

    if (length(x)==1 && abs(x-x2)<1e-10) % right boundary
        v = -ones(length(x)*length(y),1);
    elseif (length(x)==1 && abs(x-x1)<1e-10) % left boundary
        v = ones(length(x)*length(y),1);
    end
    
end