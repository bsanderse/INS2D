function u = LDC_unsteady_uBC(x,y,t,options)
% boundary conditions for u for LDC

    y2 = options.grid.y2;

    u  = zeros(length(x)*length(y),1);


    if (length(y)==1 && abs(y-y2)<1e-10)
        u = ones(length(x)*length(y),1);
    end

end