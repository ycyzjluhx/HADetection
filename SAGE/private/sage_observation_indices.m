function idx = sage_observation_indices(sys,N)
validateattributes(sys.fs,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(sys.I,{'numeric'},{'scalar','integer','positive'});
validateattributes(sys.Ta,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(sys.Tf,{'numeric'},{'scalar','real','finite','positive'});
lena = round(sys.Ta*sys.fs);
lenf = round(sys.Tf*sys.fs);
assert(lena >= 1 && lenf >= 1,'SAGE:WindowLength','Windows must contain samples.');
assert(sys.I == 1 || lenf >= lena,'SAGE:OverlappingWindows', ...
    'Overlapping windows double-count samples; use nonoverlapping windows.');
assert((sys.I-1)*lenf+lena <= N,'SAGE:WindowRange', ...
    'Observation windows exceed the available record.');
idx = reshape((1:lena)' + (0:sys.I-1)*lenf,1,[]);
end
