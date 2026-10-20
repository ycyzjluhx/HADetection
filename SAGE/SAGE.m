function [theta_out, history] = SAGE(y, u, theta, sys, options)
%SAGE Sequential profiled least-squares SAGE for complex baseband samples.
% y N-by-M; u known N-sample pilot. See README_SAGE.md for units/options.
% The original four-argument, single-output interface remains valid.
% history includes iteration zero, all parameters, and observed residual SSE.
if nargin < 5, options = struct(); end
cfg = sage_config(u, sys, options);
sage_validate_y(y, cfg);
sage_validate_theta(theta, sys, cfg);
u = u(:);
theta_out = theta(:);
for l = 1:sys.L
    theta_out(l).tau = round(theta_out(l).tau*sys.fs)/sys.fs;
end
Y = y(cfg.idx,:).'; % transpose orientation WITHOUT conjugating the data
model = zeros(size(Y));
paths = cell(sys.L,1);
for l = 1:sys.L
    paths{l} = sage_atom(theta_out(l), u, sys, cfg);
    model = model + paths{l};
end
residual = Y-model;
history.objective = nan(cfg.options.maxIter+1,1);
history.relativeImprovement = nan(cfg.options.maxIter+1,1);
history.parameterStep = nan(cfg.options.maxIter+1,1);
history.theta = cell(cfg.options.maxIter+1,1);
history.objective(1) = sum(abs(residual(:)).^2);
history.theta{1} = theta_out;
history.converged = false;
history.stopReason = 'maxIter';
dataEnergy = max(sum(abs(Y(:)).^2), realmin);
for iteration = 1:cfg.options.maxIter
    previous = theta_out;
    for l = 1:sys.L
        % beta=1: CURRENT other paths, including updates from this sweep.
        X = residual + paths{l};
        candidate = sage_fit_path(X, theta_out(l), u, sys, cfg);
        newPath = sage_atom(candidate, u, sys, cfg);
        newResidual = X-newPath;
        % Current coordinates are always candidates, so an increase is not
        % expected mathematically. Reject numerical increases as a safeguard.
        if sum(abs(newResidual(:)).^2) <= sum(abs(residual(:)).^2)
            theta_out(l) = candidate;
            paths{l} = newPath;
            residual = newResidual;
        end
    end
    cost = sum(abs(residual(:)).^2);
    lastCost = history.objective(iteration);
    relative = abs(lastCost-cost)/max(lastCost, eps*dataEnergy);
    step = sage_parameter_step(theta_out, previous, cfg);
    history.objective(iteration+1) = cost;
    history.relativeImprovement(iteration+1) = relative;
    history.parameterStep(iteration+1) = step;
    history.theta{iteration+1} = theta_out;
    if iteration >= cfg.options.minIter && ...
            relative <= cfg.options.tolObjective && step <= cfg.options.tolParameter
        history.converged = true;
        history.stopReason = 'objectiveAndParameters';
        break
    end
end
history.iterations = iteration;
history.objective = history.objective(1:iteration+1);
history.relativeImprovement = history.relativeImprovement(1:iteration+1);
history.parameterStep = history.parameterStep(1:iteration+1);
history.theta = history.theta(1:iteration+1);
history.observationIndices = cfg.idx;
history.options = cfg.options;
end
