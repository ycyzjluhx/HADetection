function [out, energy] = zfun(xel, c, u, thetal, sys)
%ZFUN z=<unit-gain path,xel> and exact observed unit-gain energy.
% Gain ML=z/energy; profiled score=abs(z)^2/energy. xel M-by-N; c M-by-1.
u = u(:); c = c(:);
assert(isequal(size(xel),[sys.M,numel(u)]),'SAGE:ObservationSize', ...
    'xel must be M-by-N.');
assert(numel(c)==sys.M,'SAGE:SteeringSize','c must contain M elements.');
idx = sage_observation_indices(sys,numel(u));
ud = sage_delayed(u,thetal.tau,sys.fs);
b = ud(idx).'.*exp(1i*2*pi*thetal.fdopp*(idx-1)/sys.fs);
out = sum(conj(b).*(c'*xel(:,idx)));
energy = sum(abs(c).^2)*sum(abs(b).^2);
end
