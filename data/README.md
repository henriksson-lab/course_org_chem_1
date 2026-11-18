# Teaching data and dictionary

All files are constructed teaching examples. None contains laboratory
measurements or student records. `python3 tools/build_data.py` regenerates them
deterministically. The table values and their meanings are part of the exercise,
not experimental reference values for identifying an unknown real sample.

| File | Columns and units | Interpretation |
| --- | --- | --- |
| `amounts.csv` | `reagent`: name; `role`: reactant; `mass_g`: grams; `molar_mass_g_mol`: g/mol; `stoich`: balanced-equation coefficient | Simulated preparation table for salicylic acid + acetic anhydride → aspirin + ethanoic acid; all coefficients 1. `mass/M` gives mol, then divide by the coefficient to compare reaction extents |
| `isolation.csv` | `sample_id`: identifier; `product`: name; `isolated_mass_g`: grams; `product_molar_mass_g_mol`: g/mol; `product_stoich`: coefficient | P1 is a simulated isolation from the amounts table; purity is unknown |
| `conformers.csv` | `conformer`: state label; `delta_G_J_mol`: relative free energy, J/mol | Simplified two-state teaching model at 298.15 K; `gauche_combined` combines the gauche ensemble and its degeneracy is already included in the supplied free energy. Not a measured conformer free energy |
| `enantiomers.csv` | `sample_id`: identifier; `R_percent`, `S_percent`: percentages | Constructed two-enantiomer mixtures; percentages sum to 100; no additional stereoisomers or impurities assumed |
| `tlc.csv` | `sample_id`: lane identity; `spot_distance_cm`, `front_distance_cm`: distances from origin in cm | One row per spot, not per independent experiment. Crude has two rows because it has two spots. These lanes share an assumed plate and eluent |
| `pka.csv` | `compound`: name; `pKa`: dimensionless approximate value; `solvent`: water; `status`: provenance label | Approximate teaching constants, not new measurements; temperature/ionic-strength differences are ignored in the model |
| `ir.csv` | `sample_id`: molecule label; `wavenumber_cm1`: cm⁻¹; `transmittance_percent`: percent | Four synthetic schematic IR traces from sums of Gaussian dips. Bands are selected to teach functional groups; shapes/intensities are illustrative, not a validated instrument simulation |
| `nmr_peaks.csv` | `sample_id`: unknown identifier; `shift_ppm`: chemical shift in ppm; `integration`: relative proton count; `multiplicity`: word label; `assignment`: answer key label | Simulated simple first-order proton NMR peak list. U1 is ethyl ethanoate (C₄H₈O₂). Hide the assignment column while solving. Shifts are illustrative; this is a peak table, not a digitised spectrum |

There are no intentional missing values. Stop and inspect the source if an
import creates missing numbers. Decimal marks are dots, fields are separated by
commas, and the first row is the header. Preserve sample labels during plotting.

Keep raw tables unchanged. Save calculations in the Rmd and derived outputs in
`output/`. In a real experiment record instrument/solvent conditions, sample
history, measurement uncertainty and any excluded observations as well.
