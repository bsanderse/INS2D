function [T] = getFOM_Temperature(RT,t,options)
% get FOM velocity based on ROM coefficients

BT   = options.rom.BT;
%Vbc = options.rom.Vbc;

% V   = B*R + Vbc;
T   = BT*RT;
