# QPSO Paper Reproduction Scripts

This folder now contains MATLAB code for the benchmark workflow in:

Alvarez-Alvarado et al., "Three novel quantum-inspired swarm optimization algorithms using different bounded potential fields", Scientific Reports 11, 11655 (2021).

## Main commands

Run the full paper-sized experiment:

```matlab
[summary, raw] = run_paper_experiments();
```

Run a quick smoke test:

```matlab
[summary, raw] = run_paper_experiments( ...
    'Functions', {'f1','f10'}, ...
    'Algorithms', {'QPSO-LR','QPSO-RM','QPSO-CS'}, ...
    'NumRuns', 3, ...
    'MaxIt', 50, ...
    'MakePlots', false);
```

Outputs are written to `results/` by default:

- `summary_all.csv`
- grouped CSV summaries by function type
- `raw_results.mat`
- convergence plots when `MakePlots` is true

The full default run uses 24 benchmark functions, 6 algorithms, 30 runs, 50 particles/agents, and 1000 iterations, so it can take a while.
