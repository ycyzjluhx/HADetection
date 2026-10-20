function [u, truth, sys, options, y] = sage_validation_fixture()
% Deterministic known channel. No Communications/Statistics toolboxes.
saved = rng; restore = onCleanup(@() rng(saved)); %#ok<NASGU>
rng(20261021,'twister');
sys = struct('M',6,'L',3,'fs',1e6,'d',0.5,'I',20,'Ta',64e-6,'Tf',1e-3);
N = (sys.I-1)*round(sys.Tf*sys.fs)+round(sys.Ta*sys.fs);
u = exp(1i*(pi/4+pi/2*randi([0 3],N,1)));
truth = struct('tau',num2cell([2 6 11]/sys.fs), ...
    'phi',num2cell([36 84 142]*pi/180), ...
    'fdopp',num2cell([-180 70 240]), ...
    'amp',num2cell([1 0.7 0.5].*exp(1i*[0.4 -1.1 2.2])));
truth = truth(:);
options = struct('tauGrid',(0:12)/sys.fs,'phiGrid',(10:2:170)*pi/180, ...
    'fdGrid',-300:10:300,'maxIter',40,'minIter',2, ...
    'tolObjective',1e-9,'tolParameter',1e-7);
y = sage_reference_signal(u,truth,sys);
end
