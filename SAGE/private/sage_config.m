function cfg = sage_config(u,sys,options)
validateattributes(u,{'double'},{'vector','nonempty','finite'});
validateattributes(sys.M,{'numeric'},{'scalar','integer','positive'});
validateattributes(sys.L,{'numeric'},{'scalar','integer','positive'});
validateattributes(sys.d,{'numeric'},{'scalar','real','finite','positive'});
cfg.N = numel(u); cfg.M = sys.M;
cfg.idx = sage_observation_indices(sys,cfg.N);
cfg.t = (cfg.idx-1)/sys.fs; % absolute times, including gaps
optionsDefault = struct('tauGrid',(0:min(floor(4e-6*sys.fs),cfg.N-1))/sys.fs, ...
    'phiGrid',linspace(0,pi,181),'fdGrid',-200:5:200, ...
    'maxIter',50,'minIter',2,'tolObjective',1e-8,'tolParameter',1e-5);
assert(isstruct(options) && isscalar(options),'SAGE:Options','options must be a scalar struct.');
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(optionsDefault,names{k}),'SAGE:OptionName','Unknown option: %s',names{k});
    optionsDefault.(names{k}) = options.(names{k});
end
options = optionsDefault;
for name = {'tauGrid','phiGrid','fdGrid'}
    key = name{1}; value = options.(key);
    validateattributes(value,{'numeric'},{'vector','real','finite','nonempty'});
    options.(key) = unique(value(:).');
end
assert(all(options.tauGrid>=0) && all(round(options.tauGrid*sys.fs)<cfg.N), ...
    'SAGE:DelayRange','Delay grid is outside the finite record.');
options.tauGrid = unique(round(options.tauGrid*sys.fs))/sys.fs;
assert(all(options.phiGrid>=0 & options.phiGrid<=pi),'SAGE:AngleRange', ...
    'Use the identifiable ULA half-plane 0 <= phi <= pi.');
assert(all(abs(options.fdGrid)<sys.fs/2),'SAGE:DopplerRange', ...
    'Doppler grid must lie strictly inside the sampling Nyquist interval.');
validateattributes(options.maxIter,{'numeric'},{'scalar','integer','positive'});
validateattributes(options.minIter,{'numeric'},{'scalar','integer','positive','<=',options.maxIter});
validateattributes(options.tolObjective,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(options.tolParameter,{'numeric'},{'scalar','real','finite','nonnegative'});
cfg.options = options;
cfg.U = zeros(numel(options.tauGrid),numel(cfg.idx));
for k = 1:numel(options.tauGrid)
    ud = sage_delayed(u,options.tauGrid(k),sys.fs);
    cfg.U(k,:) = ud(cfg.idx).';
end
assert(any(sum(abs(cfg.U).^2,2)>0),'SAGE:ZeroPilotEnergy', ...
    'No candidate has nonzero pilot energy in the observation windows.');
cfg.C = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(options.phiGrid));
cfg.D = exp(1i*2*pi*options.fdGrid(:)*cfg.t);
cfg.scale = [1/sys.fs, grid_scale(options.phiGrid),grid_scale(options.fdGrid)];
end
function scale = grid_scale(grid)
if numel(grid)>1, scale = min(diff(grid)); else, scale = 1; end
end
