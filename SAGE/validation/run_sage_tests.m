function report = run_sage_tests()
%RUN_SAGE_TESTS Regression checks of the model and confirmed implementation bugs.
addpath(fileparts(fileparts(mfilename('fullpath'))));
[u,truth,sys,opt,Y] = sage_validation_fixture();
saved = rng; restore = onCleanup(@() rng(saved)); %#ok<NASGU>
rng(2718,'twister');
names = {};
% Independent oracle detects orientation, gain phase, Doppler clock and signs.
S = zeros(sys.M,numel(u));
for l=1:sys.L, S = S+s(truth(l),u,sys); end
assert(norm(S.'-Y,'fro')/norm(Y,'fro') < 1e-12);
assert(isequal(s(truth(1),u,sys),s(truth(1),u.',sys)));
CSI = struct('delay',[truth.tau],'amp',[truth.amp], ...
    'DoA',[truth.phi],'fd',[truth.fdopp]);
[G,~] = signalGenerator(u,CSI,sys);
assert(norm(G-Y,'fro')/norm(Y,'fro') < 1e-12);
names{end+1} = 'Independent complex waveform, row pilot and Doppler clock';

X = xel(Y,2,truth,u,sys);
assert(norm(X-s(truth(2),u,sys),'fro')/norm(Y,'fro') < 1e-12);
names{end+1} = 'Hidden-data residual uses nonconjugate orientation';

idx = reshape((1:round(sys.Ta*sys.fs))'+(0:sys.I-1)*round(sys.Tf*sys.fs),1,[]);
p = truth(2); p.amp=1;
A = sage_reference_signal(u,p,sys); A = A(idx,:).';
R = randn(size(Y))+1i*randn(size(Y)); X=R(idx,:).';
c = exp(1i*2*pi*(0:sys.M-1)'*sys.d*cos(p.phi));
[z,E] = zfun(R.',c,u,p,sys);
assert(abs(z-A(:)'*X(:)) < 1e-10*max(norm(A(:))*norm(X(:)),1));
assert(abs(E-sum(abs(A(:)).^2)) < 1e-10*E);
names{end+1} = 'Hermitian matched filter and exact finite-window energy';

% Unequal complex pilot magnitudes and truncation expose global Pu errors.
one=sys; one.L=1; one.I=1; one.Ta=32/sys.fs; one.Tf=one.Ta;
v=(1:32)'.*exp(1i*(1:32)');
p=truth(1); p.tau=20/sys.fs;
o=opt; o.tauGrid=p.tau; o.phiGrid=p.phi; o.fdGrid=p.fdopp;
V=sage_reference_signal(v,p,one); seed=p; seed.amp=0;
est=SAGE(V,v,seed,one,o);
assert(abs(est.amp-p.amp)<1e-10);
names{end+1} = 'Complex gain with truncated nonconstant pilot';

seed=SAGEinit(Y,u,sys,opt);
[est,h]=SAGE(Y,u,seed,sys,opt);
[~,e]=sage_match_paths(est,truth,[1/sys.fs,2*pi/180,10,1]);
assert(max(abs(e(:,1)))<1e-12 && max(abs(e(:,2)))<1e-12);
assert(max(abs(e(:,3)))<1e-10 && max(e(:,4))<1e-7);
assert(h.objective(end)/sum(abs(Y(idx,:)).^2,'all')<1e-14);
assert(h.converged && h.iterations<opt.maxIter);
names{end+1} = 'Noiseless three-path recovery and actual early termination';

seed=truth; seed(2).fdopp=50; seed(3).fdopp=220;
[est,h2]=SAGE(Y,u,seed,sys,opt);
[~,e]=sage_match_paths(est,truth,[1/sys.fs,2*pi/180,10,1]);
assert(max(abs(e(:,3)))<1e-10 && max(e(:,4))<1e-7);
assert(h2.parameterStep(2)>0);
names{end+1} = 'Second/third path Doppler updates and sequential cancellation';

% Check that the distance reflects a changing first path, not only the last.
seed=truth; seed(1).amp=0;
o=opt; o.maxIter=1; o.minIter=1;
[~,hs]=SAGE(Y,u,seed,sys,o);
expected=max([abs([hs.theta{2}.tau]-[hs.theta{1}.tau])*sys.fs, ...
    abs([hs.theta{2}.phi]-[hs.theta{1}.phi])/(2*pi/180), ...
    abs([hs.theta{2}.fdopp]-[hs.theta{1}.fdopp])/10, ...
    abs([hs.theta{2}.amp]-[hs.theta{1}.amp])/max(norm([hs.theta{1}.amp]),1)]);
assert(abs(hs.parameterStep(2)-expected)<1e-10 && expected>0.1);
names{end+1} = 'Convergence distance includes every path and parameter';

noise=0.1*(randn(size(Y))+1i*randn(size(Y)));
[~,hn]=SAGE(Y+noise,u,SAGEinit(Y+noise,u,sys,opt),sys,opt);
assert(all(diff(hn.objective)<=1e-10*max(hn.objective(1),1)));
o=opt; o.maxIter=1; o.minIter=1; o.tolObjective=0; o.tolParameter=0;
[~,hc]=SAGE(Y+noise,u,seed,sys,o);
assert(hc.iterations==1 && ~hc.converged && strcmp(hc.stopReason,'maxIter'));
names{end+1} = 'Noisy objective monotonicity and iteration cap';

p=truth(1); p.phi=p.phi+0.003; p.fdopp=p.fdopp+1.7;
one=sys; one.L=1; V=sage_reference_signal(u,p,one);
[est,ho]=SAGE(V,u,p,one,opt);
assert(abs(est.phi-p.phi)<1e-14 && abs(est.fdopp-p.fdopp)<1e-12);
assert(ho.objective(end)<=ho.objective(1)+1e-20);
names{end+1} = 'Current off-grid candidate is retained';

bad=sys; bad.Ta=2*bad.Tf;
expect_error(@() SAGEinit(Y,u,bad,opt),'SAGE:OverlappingWindows');
bad=sys; bad.I=sys.I+1;
expect_error(@() SAGEinit(Y,u,bad,opt),'SAGE:WindowRange');
expect_error(@() SAGEinit(Y.',u,sys,opt),'SAGE:ObservationSize');
expect_error(@() SAGEinit(Y,zeros(size(u)),sys,opt),'SAGE:ZeroPilotEnergy');
expect_error(@() SAGE(Y,u,truth(1),sys,opt),'SAGE:PathCount');
names{end+1} = 'Window, dimension, zero-energy and path-count checks';

report=struct('passed',true,'count',numel(names),'names',{names}, ...
    'matlabVersion',version,'release',version('-release'));
fprintf('PASS: %d regression groups, MATLAB %s (%s)\n',report.count,version,version('-release'));
for k=1:numel(names), fprintf('  PASS %s\n',names{k}); end
end
function expect_error(f,id)
try
    f();
catch err
    assert(strcmp(err.identifier,id),'Unexpected error: %s',err.identifier);
    return
end
error('SAGE:MissingError','Expected error %s was not raised.',id);
end
