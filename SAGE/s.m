function out = s(theta, u, sys)
%S Noise-free single path, M-by-N. Angle phi is measured from ULA endfire.
% d spacing/wavelength; tau seconds; fdopp Hz; amp complex voltage gain.
u = u(:);
ud = sage_delayed(u,theta.tau,sys.fs);
t = (0:numel(u)-1)/sys.fs;
c = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(theta.phi));
out = theta.amp*c*(ud.'.*exp(1i*2*pi*theta.fdopp*t));
end
