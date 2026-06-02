function [gbest, fgbest, history] = GA_Optimizer(objfun, n, d, MaxIt, lb, ub, options)
%GA_OPTIMIZER Real-coded genetic algorithm baseline without toolbox dependency.

if nargin < 8 || isempty(options)
    options = struct();
end

lb = expand_bound(lb, d);
ub = expand_bound(ub, d);
pc = get_option(options, 'CrossoverProbability', 0.8);
pm = get_option(options, 'MutationProbability', min(0.2, 1 ./ d));
eliteCount = get_option(options, 'EliteCount', 1);
tournamentSize = get_option(options, 'TournamentSize', 2);
mutationScale = get_option(options, 'MutationScale', 0.1);

range = ub - lb;
population = repmat(lb, n, 1) + rand(n, d) .* repmat(range, n, 1);
fitness = evaluate_population(objfun, population);
[fgbest, idx] = min(fitness);
gbest = population(idx, :);
history = zeros(MaxIt, 1);

for t = 1:MaxIt
    [fitness, order] = sort(fitness);
    population = population(order, :);
    nextPopulation = zeros(n, d);
    nextPopulation(1:eliteCount, :) = population(1:eliteCount, :);

    fillIndex = eliteCount + 1;
    while fillIndex <= n
        parent1 = population(tournament_select(fitness, tournamentSize), :);
        parent2 = population(tournament_select(fitness, tournamentSize), :);

        if rand() < pc
            blend = rand(1, d);
            child1 = blend .* parent1 + (1 - blend) .* parent2;
            child2 = blend .* parent2 + (1 - blend) .* parent1;
        else
            child1 = parent1;
            child2 = parent2;
        end

        sigma = mutationScale .* (1 - (t - 1) ./ max(MaxIt - 1, 1)) .* range;
        child1 = mutate_child(child1, pm, sigma);
        child2 = mutate_child(child2, pm, sigma);
        child1 = enforce_bounds(child1, lb, ub);
        child2 = enforce_bounds(child2, lb, ub);

        nextPopulation(fillIndex, :) = child1;
        fillIndex = fillIndex + 1;
        if fillIndex <= n
            nextPopulation(fillIndex, :) = child2;
            fillIndex = fillIndex + 1;
        end
    end

    population = nextPopulation;
    fitness = evaluate_population(objfun, population);

    [candidateBest, idx] = min(fitness);
    if candidateBest < fgbest
        fgbest = candidateBest;
        gbest = population(idx, :);
    end
    history(t) = fgbest;
end
end

function idx = tournament_select(fitness, tournamentSize)
n = numel(fitness);
candidates = randi(n, tournamentSize, 1);
[~, bestLocal] = min(fitness(candidates));
idx = candidates(bestLocal);
end

function child = mutate_child(child, mutationProbability, sigma)
mask = rand(size(child)) < mutationProbability;
child(mask) = child(mask) + sigma(mask) .* randn(1, nnz(mask));
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
    error('GA_Optimizer:BadBounds', 'Bounds must be scalar or length %d.', d);
end
end

function value = get_option(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
