# Growth cost of protein over- and underexpression in nonlinear resource allocation models

## Overview

This repository contains R scripts and model files for analyzing the growth effects of protein over- and underexpression using Growth Balance Analysis (GBA). The framework explicitly accounts for nonlinear enzyme kinetics, molecular crowding, and cellular resource constraints.

We study growth costs associated with:
- Functional proteins (metabolic enzymes and molecular machines such as DNAP, RNAP, and ribosomes)
- Idle proteins without direct catalytic function
- Proteins that produce inhibitory metabolites
- Alternative reactions

We further analyze how these costs depend on environmental conditions, including external substrate concentrations and the availability of alternative metabolic routes.

The used R packages are specified in the file `renv.lock`. To recreate the environment, run `renv::restore()`.


## Repository structure

```text
├── code/
│ └── Models/
├── data/
│ └── proteomics/
├── figures/
├── README.md
├── LICENCE
└── renv.lock
```

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
  Vary the proteome fraction (`phi`) of individual proteins and predict growth rate changes.

- `test_fuel_flux_tradeoff.R`  
  Analyze tradeoffs between toxic metabolite production and efflux pump expression.

- `test_protein_cost_dekel.R`  
  Analyze the cost of LAC protein expression across environmental conditions.

- `test_protein_cost_idle_Q.R`  
  Compute growth costs associated with a single idle protein.

- `test_protein_cost_trans.R`  
  Evaluate the impact of adding an alternative transporter on the burden of an existing transporter.

- `test_ribosome_inhibition.R`  
  Reproduce ribosome inhibition growth laws following the framework of Terry Hwa.


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
- `Readmodels.R`
- `solver_loop.R` — iteration over different initial conditions
- `uni_color.R` — color definitions

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
- `M10_fuel_efflux.ods` — toxic compound (F) production with efflux pump
- `B.ods` — minimal model with three reactions


### Data

#### `data/`
Outputs of simulations.
Each file name contains the corresponding model name.

#### `data/proteomics/`
Proteomics data from *Escherichia coli* from:

- https://doi.org/10.15252/msb.20209536  
- https://doi.org/10.1038/nbt.3418  

Processed versions of these datasets are generated using the scripts in `process_proteomics.R`.

- `ecocyc_groups.txt` - GO annotation from Ecocyc


## Citation
