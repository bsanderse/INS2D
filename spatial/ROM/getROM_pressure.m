function [q] = getROM_pressure(p,t,options)
% get ROM velocity coefficients

Bp   = options.rom.Bp;
q    = Bp'*p;

