function theta_out = SAGEinit(y, u, sys, options)
%SAGEINIT Successive cancellation with joint grid maximum-likelihood seeds.
% Same finite-window complex model and grids as SAGE. Search one delay at a
% time without storing a three-dimensional dictionary of full signals.
if nargin < 4, options = struct(); end
cfg = sage_config(u, sys, options);
sage_validate_y(y, cfg);
residual = y(cfg.idx,:).';
theta_out = repmat(struct('tau',0,'phi',0,'fdopp',0,'amp',0),sys.L,1);
for l = 1:sys.L
    bestScore = -inf;
    for k = 1:numel(cfg.options.tauGrid)
        energy = sys.M*sum(abs(cfg.U(k,:)).^2);
        if energy == 0, continue; end
        z = cfg.C' * ((residual.*conj(cfg.U(k,:))) * cfg.D');
        [score, index] = max(abs(z(:)).^2/energy);
        if score > bestScore
            [ip, id] = ind2sub(size(z),index);
            bestScore = score;
            theta_out(l) = struct('tau',cfg.options.tauGrid(k), ...
                'phi',cfg.options.phiGrid(ip),'fdopp',cfg.options.fdGrid(id), ...
                'amp',z(index)/energy);
        end
    end
    residual = residual - sage_atom(theta_out(l),u(:),sys,cfg);
end
end
