function [V] = getFOM_velocity(R,t,options)
% get FOM velocity based on ROM coefficients

B   = options.rom.B;
Vbc = options.rom.Vbc;

V   = B*R + Vbc;

