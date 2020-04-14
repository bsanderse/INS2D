function [p] = getFOM_pressure(q,t,options)
% get FOM pressure based on ROM coefficients

Bp  = options.rom.Bp;
p   = Bp*q;

