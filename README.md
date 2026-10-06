# Data and driver scripts for fixed-point acceleration in polymer SCFT

This repository accompanies the manuscript **“Accelerating Polymer Self-Consistent Field Theory with Robust Preconditioned Acceleration Algorithms.”** It contains the processed numerical data, representative Julia driver scripts, and final figures used in the submitted work.

## Contents

- `data/benchmarks/`: BCC and Frank--Kasper sigma convergence traces, selected parameters, threshold crossings, parameter grids, wall-clock comparisons, and discretization results.
- `data/gyroid/`: neighboring-solution results and the mean evaluation-count parameter grid for the Gyroid tests.
- `data/hex/`: target-HEX convergence fractions, parameter-sensitivity results, and outcome statistics for the random-initialization tests.
- `scripts/benchmark/`: principal convergence and history-length timing drivers.
- `scripts/discretization/`: contour-step and spatial-resolution drivers.
- `scripts/gyroid/`: Gyroid neighboring-solution and parameter-grid drivers.
- `scripts/hex/`: random-initialization warm-up drivers for Anderson acceleration, NGMRES, and OACCEL.
- `figures/final/`: the submitted main-text and Supporting Information figures in PDF format.
- `notebooks/`: plotting notebooks retained from the original data release.

The main data-to-figure correspondence is:

| Figure | Data |
| --- | --- |
| Main Fig. 1 | `data/benchmarks/history_cost_vs_m.csv` |
| Main Figs. 2--3; SI Figs. S8--S9 | `data/benchmarks/convergence_traces/`, `selected_parameters.csv`, and `threshold_crossings.csv` |
| Main Figs. 4--5 | `data/gyroid/convergence_summary.csv` and `parameter_grid.csv` |
| Main Figs. 6--7 | `data/hex/common_vs_selected_target_fraction.csv` and `selected_parameter_outcomes.csv` |
| SI Figs. S1--S4 | `data/benchmarks/parameter_grid_bcc.csv` and `parameter_grid_sigma.csv` |
| SI Fig. S5 | `data/benchmarks/discretization_scan.csv` |
| SI Fig. S6 | `data/gyroid/mean_neval_parameter_grid.csv` |
| SI Fig. S7 | `data/hex/parameter_scan_warmup70.csv` |

## Driver scripts

The Julia files record the model construction, algorithms, parameter grids, stopping criteria, and solve calls used for the calculations. Run them from an environment containing Polyorder and its companion packages, for example with `julia --project=/path/to/project ...`.

Documentation for Polyorder is available at <https://www.yxliu.group/Polyorder.jl/dev/>.
