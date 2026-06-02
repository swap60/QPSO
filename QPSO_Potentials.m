function [gbest, fgbest, history] = QPSO_Potentials(objfun, n, d, MaxIt, lb, ub, variant, options)
%QPSO_POTENTIALS Quantum-inspired PSO using bounded potential fields.
%   VARIANT can be 'LR', 'RM', 'CS', 'QPSO-LR', 'QPSO-RM', or 'QPSO-CS'.
%
%   The particle displacement laws follow Eqs. (16), (21), and (26) from
%   Alvarez-Alvarado et al., Scientific Reports 11, 11655 (2021):
%     LR: sqrt((1-u)/u)
%     RM: sech^{-1}(sqrt(u)) = acosh(1/sqrt(u))
%     CS: log(1/u)^(2/3)

if nargin < 8 || isempty(options)
    options = struct();
end

lb = expand_bound(lb, d);
ub = expand_bound(ub, d);
variant = normalize_variant(variant);
c = get_option(options, 'AccelerationCoefficients', [2 2]);
tolerance = get_option(options, 'Tolerance', 0);

x = initialize_population(n, d, lb, ub);
pbest = x;
fpbest = evaluate_population(objfun, pbest);
[fgbest, idx] = min(fpbest);
gbest = pbest(idx, :);
history = zeros(MaxIt, 1);

for t = 1:MaxIt
    mbest = mean(pbest, 1);

    u1 = rand(n, d);
    u2 = rand(n, d);
    weight1 = c(1) .* u1;
    weight2 = c(2) .* u2;
    lambda = weight1 ./ (weight1 + weight2 + eps);
    localAttractor = lambda .* pbest + (1 - lambda) .* repmat(gbest, n, 1);

    u = max(rand(n, d), realmin);
    fV = potential_displacement(u, variant);
    direction = 2 .* (rand(n, d) > 0.5) - 1;
    relativeWidth = abs(x - repmat(mbest, n, 1));

    x = localAttractor + direction .* relativeWidth .* fV;
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

    if tolerance > 0
        convergenceGap = norm(sum(pbest, 1) - n .* gbest);
        if convergenceGap <= tolerance
            history(t+1:end) = fgbest;
            break;
        end
    end
end
end

function fV = potential_displacement(u, variant)
switch variant
    case 'LR'
        fV = sqrt((1 - u) ./ u);
    case 'RM'
        fV = acosh(1 ./ sqrt(u));
    case 'CS'
        fV = log(1 ./ u) .^ (2 ./ 3);
    otherwise
        error('QPSO_Potentials:UnknownVariant', 'Unknown QPSO variant.');
end
end

function variant = normalize_variant(variant)
variant = upper(strrep(char(variant), 'QPSO-', ''));
if ~ismember(variant, {'LR', 'RM', 'CS'})
    error('QPSO_Potentials:UnknownVariant', ...
        'Unknown variant "%s". Use LR, RM, or CS.', variant);
end
end

function x = initialize_population(n, d, lb, ub)
x = repmat(lb, n, 1) + rand(n, d) .* repmat(ub - lb, n, 1);
end

function f = evaluate_population(objfun, x)
n = size(x, 1);
f = zeros(n, 1);
for i = 1:n
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
    error('QPSO_Potentials:BadBounds', ...
        'Bounds must be scalar or length %d.', d);
end
end

function value = get_option(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
