# Camp_Century_complex_26Al10Be_modelling
Camp Century complex 26Al/10Be modelling

This MatLab script models the burial and paleo-exposure history of the Camp Century subglacial sediment given the luminescence-constrained depositional age (406 +- 35 ka) of the upper most sediment (1059-4)
during MIS 11. Use the script "lumin_cosmo_burial_model_jun22.m" to run the model.

Using the observed 26Al and 10Be concentrations and uncertainties in the upper-most (1059-4) and the lower-most (1063-7) sediments, this script first performs a Monte Carlo simulation of nuclide concentrations and
then calculates average nuclide concentrations weighted by measurement uncertainty. Then, the 26Al and 10Be concentrations of 1059-4 are corrected for 406 +- 35 kyr of burial. This yields the 26Al/10Be ratio of the upper sediment at the time of sediment deposition
in a small surface stream, assuming the sediment was sufficiently buried (by either sediment or ice) to prevent additional nuclide production since 406 +- 35 ka.The next step calculates the inventory of nuclides that could accumulate during a range of paleo-exposure periods of the upper sediment given the 26Al/10Be ratio of the upper sediment at the time of deposition. This
also constrains the maximum duration (16 kyr) of ice-free surface exposure of the upper sediment when the pre-MIS11 inherited mean 26Al/10Be ratio is equal to 0. For the lower sediment (1063-7), this script corrects the 26Al/10Be ratio for 406 kyr of burial PLUS an additional period of burial (0-16 kyr) that is equivalent to the range of exposure durations of the upper sediment, a total burial correction time of up to 422 kyr.

For both sediments, the script propagates the uncertainties of the nuclide concentrations and the luminescence age.

This script generates two plots: Figure 1 is a two-isotope plot (banana plot) of 10Be concentrations vs the 26Al/10Be ratios of the upper and lower sediments: 1) weighted mean observed values and 2)
corrected for 406 +- 35 kyr of burial. For the upper sediment, the modelled 26Al/10Be and 10Be concs given the range of MIS 11 paleo-exposure durations is shown (blues and greens). For the lower
sediment, the 26Al/10Be and 10Be concs are shown given the additional period of burial (reds and yellows).

Figure 2 contains two subplots: shows the paleo-exposure period of the uppper sediment vs the resulting 1) inherited 26Al/10Be ratio, and 2)the inherited total burial history.
