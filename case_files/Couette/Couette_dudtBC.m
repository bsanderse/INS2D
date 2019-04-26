function dudt = Couette_dudtBC(x,y,t,options)
% accelerating lid which goes to u=1

dudt = 0*ones(length(x)*length(y),1);

if (length(y)==1 && abs(y-y2)<1e-10)
    dudt = exp(-t)*ones(length(x)*length(y),1);
elseif (length(y)==1 && abs(y-y1)<1e-10)
    dudt = 0*ones(length(x)*length(y),1);
end


end