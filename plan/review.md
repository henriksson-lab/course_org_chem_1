# First teaching draft: implementation and review

19 November 2026. Scope: the new 7.5-credit Organic Chemistry 1, with a complete
five-week teaching route and a deliberately small R/Rmd component.

## Delivered

- Adapted the quantitative-biology Quarto scaffold, including sidebar, search,
  Reveal slides, speaker notes, controls and `docs/` output. The default build
  renders without staging files. Its cache stays in the writable workspace.
- Built 15 units, each with lecture, worked practice and explained self-check.
  Seven additional course/setup/reference pages bring the student source total
  to 52 pages.
- Built 240 lecture slides. Every deck has ten meaningful visual slides out of
  sixteen (62.5%). There are 150 original SVG diagrams and two original R plots.
- Supplied eight short R activities, constructed CSVs with dictionaries, a
  functional Rmd report, an RStudio project and an explicit student ZIP.
- Excluded planning, instructor files and tools from the student website and ZIP.

## Scientific review and corrections

After the scaffold and acid–base pilot, three parallel workers built separate
topic groups. Each group received an independent read-only scientific review
from another worker. Parent review included selected original SVGs and browser
screenshots. These are author/cross-unit checks, not departmental approval.

Corrections included unambiguous skeletal terminal atoms, actual lone-pair dots
in the early oxygen-charge diagram, chair-flip geometry, resonance versus
equilibrium arrows, carbocation face approaches, protonated carbonyl labels,
activation-barrier baseline, carbon NMR label spacing, and the dependency of
theoretical yield on starting amounts rather than isolated mass.

Review checked lactic-acid Fischer priorities/configuration, meso relationships,
E/Z and Newman comparisons, ring flips, elementary reaction accounting,
amide/amine and aromatic-nitrogen reasoning, and the simulated ester spectral
assignments. Teaching models explicitly label their assumptions and provenance.

## Verification

- Full Quarto rendering completes without unresolved-link warnings.
- Structural checks cover all unit pages, SVG XML, source image descriptions,
  per-deck visual quotas, rendered local targets (including Reveal lazy images),
  rendered image alternative text and public/private boundaries.
- Complete R examples from all eight relevant pages run using only the
  extracted student project; numerical anchors agree with independent answers.
- A new R process knits the report from an extracted student ZIP.
- Local Chrome checks all 240 slides with all figures preloaded: no missing
  figures or content extending beyond the slide canvas. Representative
  acid–base, stereochemistry and spectroscopy slides were visually inspected.

The machine-readable slide inventory and browser results are adjacent to this
file. Re-run `make render`, `make check` and optional `make layout` after changes.

## Classroom decisions still needed

Confirm the official 5KE231 syllabus and reading list, examination rules,
actual course sequence and teaching language. Supply approved laboratory
protocols, facilities and arrangements; current practical material supports
preparation and interpretation. Pilot pacing and projection readability with
students, especially spatial chemistry and the carbonyl block. The simulated
spectra should be complemented with real examples in classroom delivery.
