function [gbest, fgbest, history] = FFA_Optimizer(objfun, n, d, MaxIt, lb, ub, options)
%FFA_OPTIMIZER Firefly optimization baseline without toolbox dependency.

if nargin < 8 || isempty(options)
    options = struct();
end

lb = expand_bound(lb, d);
ub = expand_bound(ub, d);
alpha0 = get_option(options, 'Alpha', 0.25);
alphaDamp = get_option(options, 'AlphaDamp', 0.97);
beta0 = get_option(options, 'Beta0', 1);
gamma = get_option(options, 'Gamma', 1);

range = ub - lb;
x = repmat(lb, n, 1) + rand(n, d) .* repmat(range, n, 1);
f = evaluate_population(objfun, x);
[fgbest, idx] = min(f);
gbest = x(idx, :);
history = zeros(MaxIt, 1);

for t = 1:MaxIt
    alpha = alpha0 .* alphaDamp .^ (t - 1);
    for i = 1:n
        for j = 1:n
            if f(j) < f(i)
                scaledDiff = (x(i, :) - x(j, :)) ./ max(range, eps);
                rij2 = sum(scaledDiff .^ 2);
                beta = beta0 .* exp(-gamma .* rij2);
                randomStep = alpha .* (rand(1, d) - 0.5) .* range;
                x(i, :) = x(i, :) + beta .* (x(j, :) - x(i, :)) + randomStep;
                x(i, :) = enforce_bounds(x(i, :), lb, ub);
                f(i) = objfun(x(i, :));
            end
        end
    end

    [candidateBest, idx] = min(f);
    if candidateBest < fgbest
        fgbest = candidateBest;
        gbest = x(idx, :);
    end
    history(t) = fgbest;
end
end

function f = evaluate_population(objfun, x)
f = zeros(size(x, 1), 1);
for i = 1:size(x, 1)
    f(i) = objfun(x(i, :));
end
end

function x = enforce_bounds(x, lb, ub)
x = max(x, lb);
x = min(x, ub);
end

function bound = expand_bound(bound, d)
if isscalar(bound)
    bound = repmat(bound, 1, d);
else
    bound = bound(:)';
end
if numel(bound) ~= d
    error('FFA_Optimizer:BadBounds', 'Bounds must be scalar or length %d.', d);
end
end

function value = get_option(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
