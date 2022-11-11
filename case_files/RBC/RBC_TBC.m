function T = RBC_TBC(x,y,t,options)
% boundary conditions for u for LDC

    y1 = options.grid.y1;
    y2 = options.grid.y2;

    T  = zeros(length(x)*length(y),1);

%     % set temperature of bottom to 1
%     if (length(y)==1 && abs(y-y1)<1e-10)
%         T = ones(length(x)*length(y),1);
%     end
    
    % set temperature of bottom to 0.5
%     if (length(y)==1 && abs(y-y1)<1e-10)
%         T = 0.5*ones(length(x)*length(y),1);
%     end
    % set temperature of top to -0.5
%     if (length(y)==1 && abs(y-y2)<1e-10)
%         T = -0.5*ones(length(x)*length(y),1);
%     end
end