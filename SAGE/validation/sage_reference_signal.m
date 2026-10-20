function y = sage_reference_signal(u,truth,sys,fractionalDelay)
% Independent forward oracle, NOT s(), signalGenerator(), or private helpers.
% Optional fractional delay uses explicitly specified linear interpolation;
% it is a model-mismatch diagnostic, not a bandlimited physical channel.
if nargin<4, fractionalDelay=false; end
u = u(:); n = (0:numel(u)-1)'; y = zeros(numel(u),sys.M);
for l = 1:numel(truth)
    delay = truth(l).tau*sys.fs;
    if fractionalDelay
        delayed = interp1(n,u,n-delay,'linear',0);
    else
        delay = round(delay);
        delayed = [zeros(delay,1);u(1:end-delay)];
    end
    for m = 1:sys.M
        y(:,m) = y(:,m)+truth(l).amp*delayed.* ...
            exp(1i*2*pi*(truth(l).fdopp*n/sys.fs+(m-1)*sys.d*cos(truth(l).phi)));
    end
end
end
