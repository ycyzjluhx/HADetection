function [s, s_n] = signalGenerator( u, CSI, sys )
% Compute the noise and multipath attenuated signal as the ouput of
% channel, whose CSI is defined in 'CSI'.
%
% u: the transmitted signal.
% CSI: the predefined channel CSI including amplitude(amp), time delay(delay),
%      Doppler frequency(fd) and drection of arrival(DoA) of each path.
% sys: system information
%
% s: signal received in the receiver without noise
% s_n: noise corrupted s

L = length(CSI.delay);
M = sys.M;
d = sys.d;
fs = sys.fs;


u = u(:);
t = (0:numel(u)-1)'/fs; % do not round the duration to an integer second
s = zeros(length(u), M);
for m = 1:M
    stemp = zeros(length(u), L);
    for l = 1:L
        delayd = round(CSI.delay(l)*fs);
        ud = zeros(length(u),1);
        ud(delayd+1:end) = u(1:end-delayd);
        
        stemp(:,l) = exp(1i*2*pi*(m-1)*d*cos(CSI.DoA(l)))*CSI.amp(l)*exp(1i*2*pi*CSI.fd(l)*t).*ud;
    end
    s(:,m) = sum(stemp,2);
end

% Circular complex AWGN with E|w|^2=0.2, split between I and Q.
s_n = s + sqrt(0.2/2)*(randn(size(s))+1i*randn(size(s)));

end

