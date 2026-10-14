# Growth cost of protein over- and underexpression in nonlinear resource allocation models

bioRxiv: https://doi.org/10.64898/2026.04.21.719900

## Overview

This repository contains R scripts and model files for analyzing the growth effects of protein over- and underexpression using Growth Balance Analysis (GBA). The framework explicitly accounts for nonlinear enzyme kinetics, molecular crowding, and cellular resource constraints.

We study growth costs associated with:
- Functional proteins (metabolic enzymes and molecular machines such as DNAP, RNAP, and ribosomes)
- Idle proteins without direct catalytic function
- Proteins that produce inhibitory metabolites
- Alternative reactions

We further analyze how these costs depend on environmental conditions, including external substrate concentrations and the availability of alternative metabolic routes.

## Requirements

**R scripts.** R 4.6.1 and the packages specified in `renv.lock`. Open `prot_cost.Rproj` (the environment is activated by `.Rprofile`) and run `renv::restore()`

**Python notebook.** `toy_model_cobrak_rba.ipynb` uses a separate Python environment, defined in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate protcost
```
The kinetic (nonlinear) optimization needs the IPOPT solver, which is installed by `environment.yml`. The linear optimization in COBRA-k requires CPLEX.

## Repository structure

```text
├── code/
│ └── Models/
├── data/
│ └── experimental/
├── figures/
├── environment.yml
├── prot_cost.Rproj
├── renv.lock
├── README.md
└── LICENSE
```

## How to run

1. Open `prot_cost.Rproj`. The scripts use `here()`, so they must be run from the project root.
2. At the top of each simulation script, set
   - `modelname`: the name of a model in `code/Models/` without the file extension,
   - `is.reversible`: 1 for reversible and 0 for irreversible kinetics (`_rev` is added to the names of the output files for reversible kinetics),
   - `predict.parameters`: 0 to use the parameters in the model file.

   Some scripts loop over both kinetics.
3. Run a simulation script (`run_GBA.R` or one of the `test_*.R` scripts). The results are written to `data/`.
4. Run the corresponding plotting script. The figures are written to `figures/`. Some plotting scripts have settings at the top (for example `models_to_plot` and `relative_phi` in `plot_protein_cost_1pic.R`).

The kinetic parameters of the base model were estimated once with `predict_parameters.R`. After this step, the model file must be updated manually. The other scripts use the parameters stored in the model files.

`test_protein_cost_parallel.R` runs one reaction per core with `parallel::mclapply`. It works only on Linux and macOS, and the number of cores is set by `max_cores` and `detectCores()` in the script.

### Main analysis scripts (`code/`)

- `predict_parameters.R`  
  Predict kinetic parameters (`kcat`, `KM`). After execution, model files must be manually updated.

- `run_GBA.R`  
  Execute Growth Balance Analysis for a given model and parameter set.

- `plot_validation.R`  
  Plot ribosome-to-protein (R/P) ratios and biomass composition.

- `calculate_cost_benefit.R`  
  Quantify individual contributions to growth costs at the optimal growth state.

- `test_protein_cost_parallel.R` 
  Vary the proteome fraction (`phi`) of individual proteins and predict growth rate changes. The model is chosen by defining the `modelname` variable (must match a model in the `Models` folder without the file extension); the reversibility is defined with `is.reversible` variable (1 for reversible, 0 for irreversible).

- `test_fuel_efflux_tradeoff.R`  
  Analyze tradeoffs between toxic metabolite production and efflux pump expression.

- `test_protein_cost_dekel.R`  
  Analyze the cost of LAC protein expression across environmental conditions.

- `test_protein_cost_idle_Q.R`, `test_protein_cost_idle`  
  Compute growth costs associated with a single idle protein (with or without Q sector).

- `test_protein_cost_trans.R`  
  Evaluate the impact of adding an alternative transporter on the burden of an existing transporter.

- `test_ribosome_inhibition.R`  
  Reproduce ribosome inhibition growth laws following the framework of Terry Hwa.
  
- `toy_model_cobrak_rba.ipynb`
Compares GBA results on the toy model (`B.ods`) with RBA and with COBRA-k linear (enzyme-constrained FBA) and nonlinear optimization. It reads `data/B_rev_protein_cost.csv` and writes `data/toy_expression_scan_results_rev.json` and `figures/toy_expression_scan_rev.png`. It must be run from within the `code` folder, using the environment defined in `environment.yml`.

### Plotting scripts

Scripts for generating figures from simulation outputs:

- `plot_alt_transporter.R`
- `plot_cost_benefit_nonoptimal.R`
- `plot_dekel_simple.R`
- `plot_dekel.R`
- `plot_efflux.R`
- `plot_experimental_costs.R`
- `plot_idle.R`
- `plot_protein_cost_1pic.R`
- `plot_protein_cost_1picB.R`
- `plot_proteomics_ecocyc.R`
- `plot_ribosome_inhibition.R`
- `plot_validation.R`

### Helper scripts

Shared utility functions:

- `GBA_Exportcsv.R`
- `GBA_solver.R`
- `initialize_model.R`
- `Kinetics.R` — definitions of kinetic rate laws
- `Parameter_prediction.R` — estimation of kinetic parameters
- `process_ecocyc_annotation2.R`
- `process_proteomics.R`
- `q0_biomass.R` — linear optimization for initial solutions
- `Readmodelods.R`
- `solver_loop.R` — iteration over different initial conditions
- `uni_colors.R` — color definitions

## Model files (`code/Models/`)

The following GBA models are provided as `.ods` files:

- `M8.ods` — base model
- `M8_inh.ods` — ribosome inhibition by an external compound
- `M9_alt_trans.ods` — alternative transporter TS2
- `M9_dekel.ods` — alternative substrate utilization (lactose) by LAC proteins
- `M10_dekel_efflux.ods` — lactose utilization with efflux
- `M9_IDLE.ods` — model with an idle protein
- `M10_Q_IDLE.ods` — model with an idle protein + Q sector
- `M9_Q.ods` — constant Q sector
- `M10fuel_efflux.ods` — toxic compound (F) production with efflux pump
- `B.ods` — minimal model with three reactions

### Data

#### `data/`
Outputs of simulations. The file names contain the model name:

- `<model>_GBA.csv` — optimal growth state (`run_GBA.R`)
- `<model>_protein_cost.csv` — growth rate and cell state for forced proteome fractions (`test_protein_cost_parallel.R`)
- `<model>_tradeoff.csv` — toxic metabolite and efflux pump (`test_fuel_efflux_tradeoff.R`)
- `<model>.csv` — results of the other `test_*.R` scripts
- `_rev` in the name means reversible kinetics; without it the kinetics are irreversible
- `toy_expression_scan_results_rev.json` — results of the notebook (growth rate for forced proteome fractions with COBRA-k LP, COBRA-k NLP and RBA)

#### `data/proteomics/`
Experimental data:

- `Schmidt_proteomics_table_s6.csv`, `Schmidt_proteomics_table_s23.csv` — *Escherichia coli* proteomics from https://doi.org/10.1038/nbt.3418
- `EV2-Samples-1.csv`, `EV3-Samples-2.csv`, `EV8-AbsoluteMassFractions-1.csv`, `EV9-AbsoluteMassFractions-2.csv` — *Escherichia coli* proteomics from https://doi.org/10.15252/msb.20209536
- `ecocyc_groups.txt` — GO annotation from EcoCyc
- `scott_4B.csv` — data from Scott et al. (2010) used in `plot_idle.R`
- `biomass_calculation.xlsx` — calculation of the biomass composition

Processed versions of the proteomics datasets are generated using the scripts in `process_proteomics.R`.

## Manuscript figures

| Figure | Simulation | Plotting | Output in `figures/` |
|---|---|---|---|
| Fig. 2, Figs. S7 to S9 | `test_protein_cost_parallel.R` (`modelname = "M8"`, reversible and irreversible) | `plot_protein_cost_1pic.R` (`relative_phi = TRUE` for Figs. 2 and S7, `FALSE` for Figs. S8 and S9) | `M8_rev_rel.png`, `M8_rel.png`, `M8_rev_abs.png`, `M8_abs.png` |
| Fig. S1 | `test_protein_cost_parallel.R` (`M8`, reversible) | `plot_cost_benefit_nonoptimal.R` | `M8_rev_cost_benefit.png` |
| Table S3 | `calculate_cost_benefit.R` | | |
| Fig. 3 | `test_protein_cost_dekel.R` | `plot_dekel_simple.R` | `M9_dekel_rev_main.png` |
| Fig. S11 | `test_protein_cost_dekel.R` | `plot_dekel.R` | `M9_dekel_rev.png`, `M10_dekel_efflux_rev.png` |
| Fig. 4, Fig. S12 | `test_protein_cost_idle_Q.R`, `test_protein_cost_idle.R` | `plot_idle.R` | `M10_Q_IDLE_rev.png`, `M9_IDLE_rev.png` |
| Fig. 5 | `test_fuel_efflux_tradeoff.R` | `plot_efflux.R` | `M10fuel_efflux_rev_mu.png`, `M10fuel_efflux_rev_pump.png` |
| Fig. 6 | `test_protein_cost_trans.R` | `plot_alt_transporter.R` | `M9_alt_trans_rev.png`, `experimental_transporter.png` |
| Fig. S2 | `test_protein_cost_parallel.R` (`modelname = "B"`) | `plot_protein_cost_1pic.R` (`B` and `B_rev` in `models_to_plot`) | `B_abs.png`, `B_rev_abs.png` |
| Fig. S3 | `toy_model_cobrak_rba.ipynb` | in the notebook | `toy_expression_scan_rev.png` |
| Fig. S4 | `run_GBA.R` | `plot_validation.R` | `M8_RP_ratio.png` |
| Fig. S5 | `test_ribosome_inhibition.R` | `plot_ribosome_inhibition.R` | `M8_inh_ribosome_inhibition.png` |
| Fig. S6 | none (experimental data are in the plotting script) | `plot_experimental_costs.R` | `experimental_costs.png` |
| Fig. S10 | `test_protein_cost_parallel.R` (`modelname = "M9_Q"`) | `plot_protein_cost_1pic.R` | `M9_Q_rev_rel.png` |

## License

MIT, see `LICENSE`.

## Citation

Széliová D, Liebermeister W, Lercher M, Dourado H. Growth costs of protein allocation in models of growing cells. bioRxiv (2026). https://doi.org/10.64898/2026.04.21.719900