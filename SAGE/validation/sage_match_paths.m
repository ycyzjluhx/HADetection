function [matched, error, permutation] = sage_match_paths(estimate,truth,scales)
% Evaluation ONLY: one-to-one association minimizing scaled parameter error.
% Small-L exact permutation search avoids an assignment-toolbox dependency.
% error columns: delay [s], angle [rad], Doppler [Hz], complex gain magnitude.
L = numel(truth);
assert(numel(estimate)==L && L<=8,'SAGE:MatchingSize','Matching supports equal L <= 8.');
p = perms(1:L); scores = zeros(size(p,1),1);
for k = 1:size(p,1)
    e = errors(estimate(p(k,:)),truth);
    scores(k) = sum(sum((e./scales).^2));
end
[~,k] = min(scores); permutation = p(k,:);
matched = estimate(permutation); matched = matched(:);
error = errors(matched,truth);
end
function e = errors(a,b)
e = [[a.tau].'-[b.tau].', [a.phi].'-[b.phi].', ...
     [a.fdopp].'-[b.fdopp].', abs([a.amp].'-[b.amp].')];
end
