function u = LDC_inviscid_uBC(x,y,t,options)
% boundary conditions for u for inviscid LDC, see PhD thesis p. 70-71

    y1 = options.grid.y1;
    y2 = options.grid.y2;

    u  = zeros(length(x)*length(y),1);


    if (length(y)==1 && abs(y-y2)<1e-10)
        u = 16*(x.^2).*(1-x.^2);
    elseif (length(y)==1 && abs(y-y1)<1e-10)
        u = -ones(length(x)*length(y),1);
    end

end