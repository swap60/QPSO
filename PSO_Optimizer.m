function [gbest, fgbest, history] = PSO_Optimizer(objfun, n, d, MaxIt, lb, ub, options)
%PSO_OPTIMIZER Real-coded particle swarm optimizer baseline.

if nargin < 8 || isempty(options)
    options = struct();
end

lb = expand_bound(lb, d);
ub = expand_bound(ub, d);
c1 = get_option(options, 'CognitiveCoefficient', 2);
c2 = get_option(options, 'SocialCoefficient', 2);
wMax = get_option(options, 'InertiaMax', 0.9);
wMin = get_option(options, 'InertiaMin', 0.4);

x = repmat(lb, n, 1) + rand(n, d) .* repmat(ub - lb, n, 1);
velocityRange = ub - lb;
v = -repmat(velocityRange, n, 1) + 2 .* rand(n, d) .* repmat(velocityRange, n, 1);
vmax = 0.2 .* velocityRange;

pbest = x;
fpbest = evaluate_population(objfun, pbest);
[fgbest, idx] = min(fpbest);
gbest = pbest(idx, :);
history = zeros(MaxIt, 1);

for t = 1:MaxIt
    w = wMax - (wMax - wMin) .* (t - 1) ./ max(MaxIt - 1, 1);
    v = w .* v + c1 .* rand(n, d) .* (pbest - x) + ...
        c2 .* rand(n, d) .* (repmat(gbest, n, 1) - x);
    v = max(v, -repmat(vmax, n, 1));
    v = min(v, repmat(vmax, n, 1));

    x = x + v;
    x = enforce_bounds(x, lb, ub);

    f = evaluate_population(objfun, x);
    improved = f < fpbest;
    pbest(improved, :) = x(improved, :);
    fpbest(improved) = f(improved);

    [candidateBest, idx] = min(fpbest);
    if candidateBest < fgbest
        fgbest = candidateBest;
        gbest = pbest(idx, :);
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
x = max(x, repmat(lb, size(x, 1), 1));
x = min(x, repmat(ub, size(x, 1), 1));
end

function bound = expand_bound(bound, d)
if isscalar(bound)
    bound = repmat(bound, 1, d);
else
    bound = bound(:)';
end
if numel(bound) ~= d
    error('PSO_Optimizer:BadBounds', 'Bounds must be scalar or length %d.', d);
end
end

function value = get_option(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
