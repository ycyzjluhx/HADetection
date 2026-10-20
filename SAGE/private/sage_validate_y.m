function sage_validate_y(y,cfg)
validateattributes(y,{'double'},{'2d','finite'});
assert(isequal(size(y),[cfg.N,cfg.M]),'SAGE:ObservationSize', ...
    'y must be N-by-M, with N=numel(u).');
end
