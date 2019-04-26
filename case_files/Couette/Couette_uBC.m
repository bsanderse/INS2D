function u = Couette_uBC(x,y,t,options)
% boundary conditions for u for inviscid LDC, see PhD thesis p. 70-71

    y1 = options.grid.y1;
    y2 = options.grid.y2;

    u  = zeros(length(x)*length(y),1);


    if (length(y)==1 && abs(y-y2)<1e-10)
        % accelerating lid which goes to u=1       
        u = (1-exp(-t))*ones(length(x)*length(y),1);
    elseif (length(y)==1 && abs(y-y1)<1e-10)
        u = 0*ones(length(x)*length(y),1);
    end

end