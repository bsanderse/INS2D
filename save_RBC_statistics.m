function [Vmean,pmean,Tmean,V_var,p_var,T_var] = save_RBC_statistics(V,p,T,Vmean,pmean,Tmean,V_var,p_var,T_var)
% Vmean=options.savestatisticsRBC.Vmean;
Vmean=Vmean+V;
pmean = pmean+p;
Tmean = Tmean+T;
V_var = V_var+V.*V;
p_var = p_var+p.*p;
T_var = T_var+T.*T;
end
