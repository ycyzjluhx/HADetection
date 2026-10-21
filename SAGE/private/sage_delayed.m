function ud = sage_delayed(u,tau,fs)
% Finite-record causal integer-sample delay, zero padded, not circular.
validateattributes(tau,{'numeric'},{'scalar','real','finite','nonnegative'});
u = u(:);
delay = round(tau*fs);
assert(delay < numel(u),'SAGE:DelayRange','Delay must be shorter than the record.');
ud = zeros(size(u));
ud(delay+1:end) = u(1:end-delay);
end
