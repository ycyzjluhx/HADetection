function results = run_sage_validation(outputDir,nTrials)
%RUN_SAGE_VALIDATION Known-parameter experiments using only base MATLAB.
% From repository root:
% addpath('SAGE','SAGE/validation'); results=run_sage_validation;
% Default: 20 independent noise trials/SNR, one fixed seeded QPSK pilot.
% Outputs MAT, CSV, PNG/PDF figures and a run log identifying MATLAB release.
if nargin<1, outputDir=fullfile(fileparts(mfilename('fullpath')),'results'); end
if nargin<2, nTrials=20; end
validateattributes(nTrials,{'numeric'},{'scalar','integer','positive'});
if ~exist(outputDir,'dir'), mkdir(outputDir); end
addpath(fileparts(fileparts(mfilename('fullpath'))));
saved = rng; restore = onCleanup(@() rng(saved)); %#ok<NASGU>
diary(fullfile(outputDir,'validation_log.txt'));
closeDiary=onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('SAGE validation | %s | MATLAB %s (%s) | %s\n', ...
    datestr(now,30),version,version('-release'),computer);
tests=run_sage_tests();
[u,truth,sys,opt,Y]=sage_validation_fixture();
idx=reshape((1:round(sys.Ta*sys.fs))'+(0:sys.I-1)*round(sys.Tf*sys.fs),1,[]);
scales=[1/sys.fs,2*pi/180,10,1];
signalPower=mean(abs(Y(idx,:)).^2,'all');
rng(314159,'twister');
noiseVariance=signalPower/10^(20/10);
Yn=Y+sqrt(noiseVariance/2)*(randn(size(Y))+1i*randn(size(Y)));
initial=SAGEinit(Yn,u,sys,opt);
[estimate,history]=SAGE(Yn,u,initial,sys,opt);
[matched,errors,permutation]=sage_match_paths(estimate,truth,scales);
T=parameter_table(truth,matched);
disp(T); writetable(T,fullfile(outputDir,'parameters_20dB.csv'));
iterations=(0:history.iterations)';
C=table(iterations,history.objective,history.relativeImprovement,history.parameterStep, ...
    'VariableNames',{'Iteration','ResidualSSE','RelativeImprovement','ParameterStep'});
writetable(C,fullfile(outputDir,'convergence_20dB.csv'));

% Use the final association throughout the history, avoiding path relabeling
% at each iteration that could artificially improve apparent convergence.
errorHistory=zeros(numel(iterations),sys.L,4);
for k=1:numel(iterations)
    a=history.theta{k}(permutation);
    errorHistory(k,:,1)=abs([a.tau]-[truth.tau])*1e6;
    errorHistory(k,:,2)=abs([a.phi]-[truth.phi])*180/pi;
    errorHistory(k,:,3)=abs([a.fdopp]-[truth.fdopp]);
    errorHistory(k,:,4)=abs([a.amp]-[truth.amp]);
end
labels={'Delay absolute error (us)','Angle absolute error (deg)', ...
    'Doppler absolute error (Hz)','Complex gain absolute error'};
f=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
for j=1:4
    subplot(2,2,j); plot(iterations,errorHistory(:,:,j),'-o','LineWidth',1.2);
    grid on; xlabel('Full SAGE sweep'); ylabel(labels{j});
    legend('Path 1','Path 2','Path 3','Location','best');
    positive_axes(iterations,errorHistory(:,:,j));
end
save_figure(f,outputDir,'parameter_error_vs_iteration');
f=figure('Visible','off','Color','w','Position',[100 100 1100 450]);
subplot(1,2,1); plot(iterations,history.objective,'-o');
grid on; xlabel('Full SAGE sweep'); ylabel('Observed residual SSE');
xlim([-0.1 history.iterations+0.1]);
subplot(1,2,2); semilogy(iterations(2:end),max(history.parameterStep(2:end),1e-12),'-o');
grid on; xlabel('Full SAGE sweep'); ylabel('Maximum normalized parameter change');
title('Display floor: 10^{-12}; CSV retains exact zeros');
xlim([0.9 history.iterations+0.1]);
save_figure(f,outputDir,'convergence');

% SNR is the coherent composite signal power on the actual observed samples.
% Every trial starts from blind SIC initialization; truth is evaluation only.
snrDb=[-5 0 10 20 30]; mse=zeros(numel(snrDb),nTrials,4);
nmse=zeros(numel(snrDb),nTrials); converged=false(numel(snrDb),nTrials);
iters=zeros(numel(snrDb),nTrials); raw=[];
for q=1:numel(snrDb)
    for trial=1:nTrials
        rng(100000+1000*q+trial,'twister');
        variance=signalPower/10^(snrDb(q)/10);
        noisy=Y+sqrt(variance/2)*(randn(size(Y))+1i*randn(size(Y)));
        seed=SAGEinit(noisy,u,sys,opt);
        [est,h]=SAGE(noisy,u,seed,sys,opt);
        [~,e]=sage_match_paths(est,truth,scales);
        mse(q,trial,:)=reshape(mean(abs(e).^2,1),1,1,4);
        fitted=sage_reference_signal(u,est,sys);
        difference=fitted(idx,:)-Y(idx,:);
        nmse(q,trial)=sum(abs(difference).^2,'all')/sum(abs(Y(idx,:)).^2,'all');
        converged(q,trial)=h.converged; iters(q,trial)=h.iterations;
        assert(all(diff(h.objective)<=1e-10*max(h.objective(1),1)));
        raw=[raw;snrDb(q)*ones(sys.L,1),trial*ones(sys.L,1),(1:sys.L)',e]; %#ok<AGROW>
    end
    fprintf('SNR %+g dB: %d/%d runs met the stopping criterion.\n', ...
        snrDb(q),sum(converged(q,:)),nTrials);
end
rmse=sqrt(reshape(mean(mse,2),numel(snrDb),4));
summary=table(snrDb',rmse(:,1)*1e6,rmse(:,2)*180/pi,rmse(:,3),rmse(:,4), ...
    mean(nmse,2),mean(iters,2),mean(converged,2), ...
    'VariableNames',{'SNR_dB','Delay_RMSE_us','Angle_RMSE_deg','Doppler_RMSE_Hz', ...
    'ComplexGain_RMSE','Reconstruction_NMSE','MeanIterations','ConvergedFraction'});
disp(summary); writetable(summary,fullfile(outputDir,'rmse_vs_snr.csv'));
writetable(array2table(raw,'VariableNames',{'SNR_dB','Trial','Path','DelayError_s', ...
    'AngleError_rad','DopplerError_Hz','ComplexGainAbsError'}),fullfile(outputDir,'trial_errors.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
unitScale=[1e6,180/pi,1,1];
for j=1:4
    subplot(2,2,j); plot(snrDb,rmse(:,j)*unitScale(j),'-o','LineWidth',1.5);
    grid on; xlabel('Observed SNR (dB)'); ylabel(strrep(labels{j},'absolute error','RMSE'));
    positive_axes(snrDb,rmse(:,j)*unitScale(j));
end
save_figure(f,outputDir,'rmse_vs_snr');

% Noiseless off-grid stress test: records discretization/model error rather
% than pretending an integer-delay estimator has fractional-delay accuracy.
offgrid=truth;
delayOffset=[0.35 -0.25 0.4]; angleOffset=[0.7 -0.7 -0.6];
dopplerOffset=[2.7 3.4 -3.4];
for l=1:sys.L
    offgrid(l).tau=offgrid(l).tau+delayOffset(l)/sys.fs;
    offgrid(l).phi=offgrid(l).phi+angleOffset(l)*pi/180;
    offgrid(l).fdopp=offgrid(l).fdopp+dopplerOffset(l);
end
% Control: delay still exactly representable; only angle/Doppler are off-grid.
gridTruth=offgrid;
for l=1:sys.L, gridTruth(l).tau=truth(l).tau; end
Ygrid=sage_reference_signal(u,gridTruth,sys);
[gridEst,gridHistory]=SAGE(Ygrid,u,SAGEinit(Ygrid,u,sys,opt),sys,opt);
[gridMatched,gridErrors]=sage_match_paths(gridEst,gridTruth,scales);
writetable(parameter_table(gridTruth,gridMatched),fullfile(outputDir,'angle_doppler_offgrid.csv'));
Yoff=sage_reference_signal(u,offgrid,sys,true);
[offEst,offHistory]=SAGE(Yoff,u,SAGEinit(Yoff,u,sys,opt),sys,opt);
[offMatched,offErrors]=sage_match_paths(offEst,offgrid,scales);
writetable(parameter_table(offgrid,offMatched),fullfile(outputDir,'offgrid_noiseless.csv'));
results=struct('tests',tests,'sys',sys,'options',opt,'truth',truth,'pilot',u, ...
    'snrDb',snrDb,'nTrials',nTrials,'estimate20dB',estimate,'errors20dB',errors, ...
    'history20dB',history,'errorHistory',errorHistory,'mse',mse,'nmse',nmse, ...
    'converged',converged,'iterations',iters,'summary',summary, ...
    'angleDopplerOffgridTruth',gridTruth,'angleDopplerOffgridEstimate',gridEst, ...
    'angleDopplerOffgridErrors',gridErrors,'angleDopplerOffgridHistory',gridHistory, ...
    'offgridTruth',offgrid,'offgridEstimate',offEst,'offgridErrors',offErrors, ...
    'offgridHistory',offHistory,'matlabVersion',version,'release',version('-release'));
save(fullfile(outputDir,'sage_validation_results.mat'),'results');
fprintf('Validation complete. Artifacts: %s\n',outputDir);
end
function T=parameter_table(truth,estimate)
T=table((1:numel(truth))',[truth.tau]'*1e6,[estimate.tau]'*1e6, ...
    [truth.phi]'*180/pi,[estimate.phi]'*180/pi,[truth.fdopp]',[estimate.fdopp]', ...
    real([truth.amp])',imag([truth.amp])',real([estimate.amp])',imag([estimate.amp])', ...
    abs([estimate.amp]'-[truth.amp]'), ...
    'VariableNames',{'Path','TrueDelay_us','EstimatedDelay_us','TrueAngle_deg', ...
    'EstimatedAngle_deg','TrueDoppler_Hz','EstimatedDoppler_Hz','TrueGainReal', ...
    'TrueGainImag','EstimatedGainReal','EstimatedGainImag','ComplexGainAbsError'});
end
function save_figure(f,folder,name)
% Explicit label margins also work with headless Linux font metrics.
drawnow;
ax=findall(f,'Type','axes');
position=vertcat(ax.Position);
[~,order]=sortrows([-position(:,2),position(:,1)],[1 2]);
ax=ax(order);
if numel(ax)==4
    slots=[0.13 0.61 0.32 0.30;0.64 0.61 0.32 0.30; ...
           0.13 0.16 0.32 0.30;0.64 0.16 0.32 0.30];
else
    slots=[0.13 0.22 0.31 0.67;0.65 0.22 0.31 0.67];
end
for k=1:numel(ax)
    ax(k).Units='normalized'; ax(k).Position=slots(k,:);
    ax(k).FontSize=10;
    ax(k).YLabel.Units='normalized'; ax(k).YLabel.Position=[-0.24 0.5 0];
    ax(k).XLabel.Units='normalized'; ax(k).XLabel.Position=[0.5 -0.20 0];
end
drawnow;
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',160);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');
close(f);
end
function positive_axes(x,y)
% Error magnitudes cannot be negative; pad endpoints so markers are visible.
top=max(y(:)); if top==0, top=1; end
ylim([0,1.15*top]);
padding=max(0.1,0.03*(max(x)-min(x)));
xlim([min(x)-padding,max(x)+padding]);
end
