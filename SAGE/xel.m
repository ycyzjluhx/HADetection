function out = xel(y, l, theta, u, sys)
%XEL beta=1 hidden-data mean; y N-by-M, out M-by-N.
assert(isequal(size(y),[numel(u),sys.M]),'SAGE:ObservationSize', ...
    'y must be N-by-M, with N=numel(u).');
assert(l >= 1 && l <= sys.L && l == fix(l),'SAGE:PathIndex','Invalid path index.');
out = y.';
for k = 1:sys.L
    if k ~= l, out = out-s(theta(k),u,sys); end
end
end
