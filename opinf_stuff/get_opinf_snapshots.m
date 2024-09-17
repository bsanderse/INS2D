function [A,A_dot] = get_opinf_snapshots(A_raw,dt)
% extend this function by adding inputs such as:
% - exact FOM time derivative snapshots
% - specification, what kind of time derivative approximation should be
% used

A_dot = (A_raw(:,2:end,:) - A_raw(:,1:end-1,:))/dt;

A = A_raw(:,1:end-1,:); % explict Euler
% A = A_raw(:,2:end,:); % implicit Euler