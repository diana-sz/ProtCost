# Code

R code for growth balance analysis (GBA) of the models in `Models/`. The scripts simulate how the expression of single proteins affects the growth rate, and produce the figures of the manuscript.

## Requirements

R packages: `here`, `readODS`, `nloptr`, `lpSolve`, `Matrix`, `parallel`, `RColorBrewer`, `scales`, `viridisLite`, `xtable`, `dplyr`.

`test_protein_cost_parallel.R` uses `mclapply()`, which runs in parallel only on Linux and macOS.

All paths are built with `here()`, relative to the project root. The project root must contain the folders `code/`, `data/` (with `data/experimental/`) and `figures/` (with `figures/source_data_figures/`).

## Workflow

1. Run a simulation script. It reads a model from `Models/` and writes a CSV file to `data/`.
2. Run the matching plotting script. It reads the CSV file and writes the figure as TIF and PNG to `figures/`. For main figures, it also writes the source data to `figures/source_data_figures/`.

Settings such as the model name and reversibility are set at the top of each script. Reversible runs add `_rev` to the model name in the output file names.

## Core files

These files are sourced by the other scripts and are not run directly.

| File | Purpose |
|---|---|
| `initialize_model.R` | Reads the model, sets solver options, finds the initial solution |
| `Readmodelods.R` | Reads the `.ods` model file |
| `Kinetics.R` | Rate law: turnover times `tau()` and their derivatives `dtau()` |
| `q0_biomass.R` | Initial solution from a linear program |
| `Parameter_prediction.R` | Predicts kinetic parameters from measured proteome fractions |
| `GBA_solver.R` | Maximizes the growth rate for each condition |
| `solver_loop.R` | Calls the solver again from other starting points if it fails |
| `GBA_Exportcsv.R` | Writes the optimal states to `data/<model>_GBA.csv` |
| `uni_colors.R` | Color definitions |

## Simulation scripts

| Script | Model | Output in `data/` |
|---|---|---|
| `run_GBA.R` | M8 | `M8_GBA.csv`, `M8_rev_GBA.csv` |
| `test_protein_cost_parallel.R` | M8, M9_Q, toy | `<model>_protein_cost.csv` |
| `test_protein_cost_idle.R` | M9_IDLE | `M9_IDLE.csv`, `M9_IDLE_rev.csv` |
| `test_protein_cost_idle_Q.R` | M10_Q_IDLE | `M10_Q_IDLE.csv`, `M10_Q_IDLE_rev.csv` |
| `test_protein_cost_dekel.R` | M9_dekel, M10_dekel_efflux | `<model>_rev.csv` |
| `test_fuel_efflux_tradeoff.R` | M10fuel_efflux | `M10fuel_efflux_rev_tradeoff.csv` |
| `test_protein_cost_trans.R` | M9_alt_trans | `M9_alt_trans.csv`, `M9_alt_trans_rev.csv` |
| `test_ribosome_inhibition.R` | M8_inh | `M8_inh.csv`, `M8_inh_rev.csv` |
| `predict_parameters.R` | M8 | Prints predicted parameters as LaTeX tables |
| `calculate_cost_benefit.R` | M8 | Prints cost and benefit terms at the optimum |

`test_protein_cost_parallel.R` and `test_protein_cost_dekel.R` run one model per call. Change the model name (and `is.reversible`) at the top of the script to produce all files.

## Plotting scripts

| Script | Figure |
|---|---|
| `plot_protein_cost_1pic.R` | Figure 2, Appendix Figures S7–S10 |
| `plot_idle.R` | Figure 3, Appendix Figure S12 |
| `plot_dekel_simple.R` | Figure 4 |
| `plot_efflux.R` | Figure 5 |
| `plot_alt_transporter.R` | Figure 6 |
| `plot_cost_benefit_nonoptimal.R` | Appendix Figure S1 |
| `plot_protein_cost_1pic_toy.R` | Appendix Figure S2 |
| `plot_validation.R` | Appendix Figure S4 |
| `plot_ribosome_inhibition.R` | Appendix Figure S5 |
| `plot_experimental_costs.R` | Appendix Figure S6 |
| `plot_dekel.R` | Appendix Figure S11 |
| `plot_proteomics_ecocyc.R` | `proteomics_ecocyc.pdf` (uses `process_proteomics.R` and `process_ecocyc_annotation2.R`) |

In `plot_protein_cost_1pic.R`, `relative_phi <- TRUE` gives Figure 2, S7 and S10, and `relative_phi <- FALSE` gives S8 and S9.

Only simulations with nloptr convergence code 4 are plotted.