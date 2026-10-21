function A = sage_atom(theta,u,sys,cfg)
ud = sage_delayed(u,theta.tau,sys.fs);
b = ud(cfg.idx).'.*exp(1i*2*pi*theta.fdopp*cfg.t);
c = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(theta.phi));
A = theta.amp*c*b;
end
