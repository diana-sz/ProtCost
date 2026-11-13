# Models

Each model is one `.ods` spreadsheet. `Readmodelods.R` reads it and `Kinetics.R` builds the rate laws from it.

| Model | Reactions | Used in |
|---|---|---|
| `toy.ods` | TS, AAS, R | Appendix Figure S2 |
| `M8.ods` | TS, ATPS, AAS, NTS, ADPS, RNAP, DNAP, R | Figure 2, Appendix Figures S1, S4, S7–S9 |
| `M8_inh.ods` | as M8, with ribosome inhibitor `x_INH` | Appendix Figure S5 |
| `M9_Q.ods` | M8 + Q (housekeeping proteins) | Appendix Figure S10 |
| `M9_IDLE.ods` | M8 + IDLE (unneeded protein) | Appendix Figure S12 |
| `M10_Q_IDLE.ods` | M8 + IDLE + Q | Figure 3 |
| `M9_dekel.ods` | M8 + LAC (lactose uptake, `x_LAC`) | Figure 4, Appendix Figure S11 |
| `M10_dekel_efflux.ods` | M9_dekel + EFFLUX | Appendix Figure S11 |
| `M10fuel_efflux.ods` | M8 + FUEL + EFFLUX | Figure 5 |
| `M9_alt_trans.ods` | M8 + TS2 (alternative transporter) | Figure 6 |

## Units

Concentrations, Michaelis constants and the cell density `rho` are in g/L. Turnover numbers are in 1/h, as mass of products per mass of protein per hour.

## General rules

- Sheet names must match exactly, including upper and lower case.
- Sheets that are not listed below are ignored.
- Reactants are external (name starts with `x_`, outside the cell) or internal (inside the cell).

## Required sheets

### `M`: mass fraction matrix

One row per reactant and one column per reaction. Cell A1 is empty, the first column holds the reactant names and the first row holds the reaction names.

- The external reactants (`x_...`) must be the first rows.
- The last row must be the total protein `P`.
- The last column must be the ribosome `R`.
- Each entry is the mass fraction of a reactant in the reaction. Substrates are negative and sum to −1; products are positive and sum to +1.
- Over all rows, including external reactants, each column sums to 0 (mass is conserved).

Example (M8, shortened):

|     | TS | ATPS  | … | R     |
|-----|----|-------|---|-------|
| x_S | -1 | 0     | … | 0     |
| x_W | 0  | 0.02  | … | 0     |
| S   | 1  | -0.17 | … | 0.16  |
| …   | …  | …     | … | …     |
| P   | 0  | 0     | … | 0.05  |

### `kcat`: turnover numbers

First column: row labels. Row 1 holds the forward turnover numbers (`kcat_f`), row 2 the backward ones (`kcat_b`). One column per reaction, in the same order as in `M` (read by position). In irreversible runs, all backward values are set to 0.

|        | TS   | ATPS | … | R     |
|--------|------|------|---|-------|
| kcat_f | 55.1 | 970  | … | 108.4 |
| kcat_b | 11   | 194  | … | 0     |

### `conditions`: growth conditions

Each column after the first is one growth condition; the header row holds the condition numbers. The first column holds the row names:

- Row 1 must be `rho`, the cell density.
- The following rows give the concentration of each external reactant. The names must match the `x_` rows of `M`; their order can differ.

|     | 1     | 2   | … |
|-----|-------|-----|---|
| rho | 186   | 186 | … |
| x_S | 10    | 8   | … |
| x_W | 1e-5  | 1e-5| … |

A concentration of exactly 0 for a substrate of a reaction gives a zero rate for that reaction, and the solver fails if that reaction is used. Use a small positive value instead, for example 1e-5.

### `biomass`: lower bounds of the biomass composition

Two columns with the headers `metabolite` and `mass_fraction`. One row per internal reactant, including `P`. The names must match `M`; their order can differ. These values are only used to compute the initial solution.

## Optional sheets

### `K`, `KA`, `KI`: kinetic constants

Same layout as `M`: same rows and columns in the same order (read by position).

- `K`: Michaelis constant of each substrate and product. A value of 0 means the reaction does not depend on that reactant. Product constants are only used in reversible runs. If the sheet is missing, every substrate gets K = 0.1 g/L.
- `KA`: activation constant. A value above 0 makes the reactant an activator of the reaction; 0 means no activation.
- `KI`: inhibition constant. A value above 0 makes the reactant an inhibitor of the reaction; 0 means no inhibition.

### `min_c`, `max_c`: bounds on internal concentrations

Header row, then one row per internal reactant in the same order as in `M` (read by position). The first column holds the names and the second column the values in g/L. Defaults: `min_c` = 0, `max_c` = `rho` of the current condition.

### `min_q`, `max_q`, `min_phi`, `max_phi`: bounds per reaction

A header row with the reaction names and one row of values, one per reaction in the same order as in `M`. There is no label column.

- `min_q`, `max_q`: bounds on the flux fractions q. Defaults 0 and 100.
- `min_phi`, `max_phi`: bounds on the proteome fractions. Defaults 0 and 1.

Example (`min_phi` in M9_Q, which fixes at least 55% of the proteome as Q):

| TS | ATPS | … | Q    | R |
|----|------|---|------|---|
| 0  | 0    | … | 0.55 | 0 |

### `phi`: measured proteome fractions

Header row `mu` followed by the reaction names; one row with the measured growth rate and the proteome fraction of each reaction. Only used when `predict.parameters > 0`.

### `Info`

Free text describing the model. It is not read by the code.