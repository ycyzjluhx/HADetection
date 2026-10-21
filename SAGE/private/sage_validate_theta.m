function sage_validate_theta(theta,sys,cfg)
assert(isstruct(theta) && numel(theta)==sys.L,'SAGE:PathCount', ...
    'theta must contain exactly sys.L paths.');
for l = 1:sys.L
    for name = {'tau','phi','fdopp','amp'}
        validateattributes(theta(l).(name{1}),{'numeric'},{'scalar','finite'});
    end
    assert(isreal(theta(l).tau) && theta(l).tau>=0 && ...
        round(theta(l).tau*sys.fs)<cfg.N,'SAGE:DelayRange','Invalid initial delay.');
    assert(isreal(theta(l).phi) && theta(l).phi>=0 && theta(l).phi<=pi, ...
        'SAGE:AngleRange','Initial angle must lie in [0,pi].');
    assert(isreal(theta(l).fdopp) && abs(theta(l).fdopp)<sys.fs/2, ...
        'SAGE:DopplerRange','Invalid initial Doppler.');
end
end
