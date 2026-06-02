function [summary, raw] = run_paper_experiments(varargin)
%RUN_PAPER_EXPERIMENTS Reproduce the benchmark workflow from the QPSO paper.
%
% Full paper-size run:
%   [summary, raw] = run_paper_experiments();
%
% Fast smoke run:
%   [summary, raw] = run_paper_experiments('Functions', {'f1','f10'}, ...
%       'Algorithms', {'QPSO-LR','QPSO-RM','QPSO-CS'}, ...
%       'NumRuns', 3, 'MaxIt', 50, 'MakePlots', false);

p = inputParser();
addParameter(p, 'Functions', 'all');
addParameter(p, 'Algorithms', {'QPSO-LR','QPSO-RM','QPSO-CS','PSO','FFO','GA'});
addParameter(p, 'NumRuns', 30);
addParameter(p, 'PopulationSize', 50);
addParameter(p, 'MaxIt', 1000);
addParameter(p, 'Seed', 2021);
addParameter(p, 'ResultsDir', fullfile(pwd, 'results'));
addParameter(p, 'MakePlots', true);
parse(p, varargin{:});
cfg = p.Results;

defs = benchmark_definitions();
functionIds = normalize_function_list(cfg.Functions, defs);
algorithms = normalize_algorithm_list(cfg.Algorithms);

if ~exist(cfg.ResultsDir, 'dir')
    mkdir(cfg.ResultsDir);
end

totalRows = numel(functionIds) .* numel(algorithms);
summaryRows = cell(totalRows, 14);
raw = struct('Function', {}, 'FunctionName', {}, 'FunctionType', {}, ...
    'Algorithm', {}, 'BestValues', {}, 'Histories', {}, ...
    'BestPosition', {}, 'RunTimes', {});
row = 0;

fprintf('Running %d function(s), %d algorithm(s), %d run(s) each.\n', ...
    numel(functionIds), numel(algorithms), cfg.NumRuns);

for fIdx = 1:numel(functionIds)
    funcId = functionIds{fIdx};
    [objfun, info] = objective_fun(funcId);
    d = info.Dimension;
    lb = expand_bound(info.LowerBound, d);
    ub = expand_bound(info.UpperBound, d);

    for aIdx = 1:numel(algorithms)
        algorithm = algorithms{aIdx};
        bestValues = zeros(cfg.NumRuns, 1);
        histories = zeros(cfg.NumRuns, cfg.MaxIt);
        runTimes = zeros(cfg.NumRuns, 1);
        bestPosition = [];
        bestAcrossRuns = inf;

        fprintf('  %s / %s ... ', funcId, algorithm);
        for r = 1:cfg.NumRuns
            if ~isempty(cfg.Seed)
                rng(cfg.Seed + fIdx * 100000 + aIdx * 1000 + r, 'twister');
            end

            tic;
            [candidatePosition, candidateBest, history] = run_optimizer( ...
                algorithm, objfun, cfg.PopulationSize, d, cfg.MaxIt, lb, ub);
            runTimes(r) = toc;

            bestValues(r) = candidateBest;
            histories(r, :) = history(:)';
            if candidateBest < bestAcrossRuns
                bestAcrossRuns = candidateBest;
                bestPosition = candidatePosition;
            end
        end

        metrics = summarize_runs(bestValues, histories, info.FMin, runTimes);
        fprintf('best %.4g, accuracy %.4g, precision %.4g\n', ...
            min(bestValues), metrics.Accuracy, metrics.Precision);

        row = row + 1;
        raw(row).Function = funcId;
        raw(row).FunctionName = info.Name;
        raw(row).FunctionType = info.Type;
        raw(row).Algorithm = algorithm;
        raw(row).BestValues = bestValues;
        raw(row).Histories = histories;
        raw(row).BestPosition = bestPosition;
        raw(row).RunTimes = runTimes;

        summaryRows(row, :) = {funcId, info.Name, info.Type, algorithm, ...
            d, info.FMin, min(bestValues), mean(bestValues), std(bestValues), ...
            metrics.Accuracy, metrics.Precision, metrics.SearchSpeed, ...
            metrics.SearchAcceleration, mean(runTimes)};
    end
end

summaryRows = summaryRows(1:row, :);
summary = cell2table(summaryRows, 'VariableNames', { ...
    'Function', 'FunctionName', 'FunctionType', 'Algorithm', ...
    'Dimension', 'FMin', 'Best', 'MeanBest', 'StdBest', ...
    'Accuracy', 'Precision', 'SearchSpeed', ...
    'SearchAcceleration', 'MeanTimeSeconds'});

write_outputs(summary, raw, cfg.ResultsDir, cfg.MakePlots);
end

function [gbest, fgbest, history] = run_optimizer(algorithm, objfun, n, d, MaxIt, lb, ub)
switch upper(algorithm)
    case {'QPSO-LR', 'LR'}
        [gbest, fgbest, history] = QPSO_Potentials(objfun, n, d, MaxIt, lb, ub, 'LR');
    case {'QPSO-RM', 'RM'}
        [gbest, fgbest, history] = QPSO_Potentials(objfun, n, d, MaxIt, lb, ub, 'RM');
    case {'QPSO-CS', 'CS'}
        [gbest, fgbest, history] = QPSO_Potentials(objfun, n, d, MaxIt, lb, ub, 'CS');
    case 'PSO'
        [gbest, fgbest, history] = PSO_Optimizer(objfun, n, d, MaxIt, lb, ub);
    case {'FFO', 'FFA', 'FIREFLY'}
        [gbest, fgbest, history] = FFA_Optimizer(objfun, n, d, MaxIt, lb, ub);
    case 'GA'
        [gbest, fgbest, history] = GA_Optimizer(objfun, n, d, MaxIt, lb, ub);
    otherwise
        error('run_paper_experiments:UnknownAlgorithm', ...
            'Unknown algorithm "%s".', algorithm);
end
end

function metrics = summarize_runs(bestValues, histories, fmin, runTimes)
meanBest = mean(bestValues);
stdBest = std(bestValues);
meanHistory = mean(histories, 1);

metrics.Accuracy = abs(fmin - meanBest);
if abs(meanBest) < eps
    if stdBest < eps
        metrics.Precision = 0;
    else
        metrics.Precision = inf;
    end
else
    metrics.Precision = abs(stdBest ./ meanBest);
end

speedSeries = diff(meanHistory);
accelerationSeries = diff(speedSeries);
metrics.SearchSpeed = sqrt(0.5 .* mean(speedSeries .^ 2));
if isempty(accelerationSeries)
    metrics.SearchAcceleration = 0;
else
    metrics.SearchAcceleration = sqrt(0.5 .* mean(accelerationSeries .^ 2));
end
metrics.MeanTimeSeconds = mean(runTimes);
end

function write_outputs(summary, raw, resultsDir, makePlots)
writetable(summary, fullfile(resultsDir, 'summary_all.csv'));

groups = unique(summary.FunctionType);
for i = 1:numel(groups)
    mask = strcmp(summary.FunctionType, groups{i});
    fileName = sprintf('summary_%s.csv', lower(groups{i}));
    writetable(summary(mask, :), fullfile(resultsDir, fileName));
end

save(fullfile(resultsDir, 'raw_results.mat'), 'raw', 'summary');

if makePlots
    plot_convergence(raw, resultsDir);
end
end

function plot_convergence(raw, resultsDir)
functions = unique({raw.Function});
for fIdx = 1:numel(functions)
    funcId = functions{fIdx};
    mask = strcmp({raw.Function}, funcId);
    rows = raw(mask);
    [~, info] = objective_fun(funcId);

    fig = figure('Visible', 'off');
    hold on;
    for r = 1:numel(rows)
        meanHistory = mean(rows(r).Histories, 1);
        errorHistory = abs(meanHistory - info.FMin) + eps;
        plot(errorHistory, 'LineWidth', 1.4, 'DisplayName', rows(r).Algorithm);
    end
    hold off;
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('Iteration');
    ylabel('Mean absolute error to fmin');
    title(sprintf('Convergence behaviour: %s', funcId));
    legend('Location', 'best');
    saveas(fig, fullfile(resultsDir, sprintf('convergence_%s.png', funcId)));
    close(fig);
end
end

function functionIds = normalize_function_list(functions, defs)
allIds = {defs.Id};
if ischar(functions)
    if strcmpi(functions, 'all')
        functionIds = allIds;
    else
        functionIds = {functions};
    end
elseif isstring(functions)
    if isscalar(functions) && strcmpi(char(functions), 'all')
        functionIds = allIds;
    else
        functionIds = cellstr(functions);
    end
else
    functionIds = functions;
end

for i = 1:numel(functionIds)
    if ~ismember(lower(functionIds{i}), lower(allIds))
        error('run_paper_experiments:UnknownFunction', ...
            'Unknown function "%s".', functionIds{i});
    end
end
end

function algorithms = normalize_algorithm_list(algorithms)
if ischar(algorithms) || isstring(algorithms)
    algorithms = cellstr(algorithms);
end
for i = 1:numel(algorithms)
    algorithms{i} = upper(char(algorithms{i}));
    if strcmp(algorithms{i}, 'FFA')
        algorithms{i} = 'FFO';
    end
end
end

function bound = expand_bound(bound, d)
if isscalar(bound)
    bound = repmat(bound, 1, d);
else
    bound = bound(:)';
end
end
