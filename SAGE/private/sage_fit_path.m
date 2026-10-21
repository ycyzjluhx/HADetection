function theta = sage_fit_path(X,theta,u,sys,cfg)
% Profile complex gain at every coordinate. Keep the current coordinate
% FIRST so ties preserve it, including valid initial off-grid coordinates.
c = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(theta.phi));
ud = sage_delayed(u,theta.tau,sys.fs);
Q = [ud(cfg.idx).';cfg.U].*exp(1i*2*pi*theta.fdopp*cfg.t);
z = conj(Q)*(c'*X).';
energy = sys.M*sum(abs(Q).^2,2);
[k,gain] = best(z,energy);
grid = [theta.tau,cfg.options.tauGrid];
theta.tau = grid(k); theta.amp = gain;

ud = sage_delayed(u,theta.tau,sys.fs);
b = ud(cfg.idx).'.*exp(1i*2*pi*theta.fdopp*cfg.t);
C = [c,cfg.C];
z = C'*(X*conj(b.'));
energy = sum(abs(C).^2,1).'*sum(abs(b).^2);
[k,gain] = best(z,energy);
grid = [theta.phi,cfg.options.phiGrid];
theta.phi = grid(k); theta.amp = gain;

c = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(theta.phi));
Q = [exp(1i*2*pi*theta.fdopp*cfg.t);cfg.D].*ud(cfg.idx).';
z = conj(Q)*(c'*X).';
energy = sys.M*sum(abs(Q).^2,2);
[k,gain] = best(z,energy);
grid = [theta.fdopp,cfg.options.fdGrid];
theta.fdopp = grid(k); theta.amp = gain;
end
function [k,gain] = best(z,energy)
score = -inf(size(energy));
valid = energy>0;
score(valid) = abs(z(valid)).^2./energy(valid);
[value,k] = max(score);
assert(isfinite(value),'SAGE:ZeroPilotEnergy','No valid unit-gain atom.');
gain = z(k)/energy(k);
end
