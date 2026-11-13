Source data for Figure 6

6A/Figure_6A.csv: simulations with GBA model M9_alt_trans_rev (reversible kinetics).
  Only converged simulations with phi_TS1_relative <= 3 are included, as in the figure.
  kcat_TS2: turnover number of the alternative transporter TS2 [1/h] (line color)
  Km_TS2: Michaelis constant of TS2 for S_ext [g/L] (line type)
  phi_TS1_relative: proteome fraction of TS1 relative to its optimum (phi/phi*)
  relative_growth_rate: growth rate relative to the optimum

6B/Figure_6B_van_der_Vlag_1994.csv: experimental data (points) from
  van der Vlag et al. (1994), redrawn from Bruggeman et al. (2023).
  PTS IIA(Glc) in Salmonella typhimurium grown on maltose.
  titrated_level: titrated level / wild-type level of the transporter
  relative_growth_rate: growth rate relative to the wild type

6B/Figure_6B_fit_Bruggeman_2023.csv: phenomenological fit (line) from
  Bruggeman et al. (2023):
  relative_growth_rate = 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x,
  with x = titrated_level.
