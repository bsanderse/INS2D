function [Sfd freq phi] = PSDfunc(x, t)
% function [Sfd freq] = PSDfunc(x, t)
% 
% inputs: signal x and time signal t
% outputs: PSD of x and frequency vector
% 
% This function returns and draws the PSD of signal x. x and t should be of
% equal size. The PSD is not given in dB. To convert the PSD to dB, do the
% following:
% Sfd_dB = 20*log10(PSD)
%
% By Mareijn Willems, july 2011

if size(x) ~=size(t)
    x = x';
    if length(x) ~= length(t)
        error('x and t should be of equal size!')
    end
end

N = length(x);
T = t(end)-t(1);
dt = t(2)-t(1);
w0 = (2*pi)/T;

Sfd = conj(dt*fft(x)).*(dt*fft(x))/N;
phi = atan2(imag(fft(x)),real(fft(x)));
phi = phi(2:N/2+1);
Sfd = Sfd(2:N/2+1);
Sfd = sqrt(4*N*Sfd/(T^2));


freq = [1:1:length(Sfd)]'*w0;

% loglog(freq,Sfd); title('Sfd'); 
% plot(freq,Sfd);
% semilogx(freq,20*log10(PSD)); title('Sfd'); 
% xlabel('frequency [rad/s]'); ylabel('power [dB]')
