Source data for Figure 2

Simulations: GBA model M8. Each row of the figure varies the proteome fraction of one reaction.
Only converged simulations with phi_relative <= 5 are included, as in the figure.
The point (0, 0) added to the growth rate lines for display is not included.

2A/Figure_2A_reversible.csv: reversible kinetics (wine lines, proteome and biomass panels)
  reaction: reaction whose proteome fraction is varied
  phi: proteome fraction of that reaction
  phi_relative: phi relative to its value at the maximal growth rate (phi/phi*)
  growth_rate: growth rate [1/h]
  growth_rate_relative: growth rate relative to its maximum (mu/mu*)
  proteome_<X>: mass fraction of protein X in the proteome (each row sums to 1)
  biomass_<X>: mass fraction of component X in the biomass (each row sums to 1)

2A/Figure_2A_irreversible.csv: irreversible kinetics (grey lines)
  Columns reaction to growth_rate_relative as above, for TS, ATPS, ADPS, NTS and AAS.
