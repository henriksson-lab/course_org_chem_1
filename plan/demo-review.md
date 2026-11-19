# Daily interactive explorations: implementation and checks

25 daily demos are embedded in Exercises and linked directly from the left
sidebar under the corresponding Week and Day, alongside the lecture. The catalogue is `demos/index.qmd`;
`demos.json` maps the days to exercise anchors, shared partials, compatibility
URLs, titles and runtime groups. The old demo URLs redirect to Exercises.

## Teaching purpose

Each page contains a prediction, a short exploration and an explanation prompt,
with all necessary definitions and model assumptions supplied locally. The
activities fit an existing 5–10 minute example/discussion block or independent
practice; they do not add a new mandatory lecture block or displace the paper
and R work.

The early demos explore representations, hydrogen counts, functional groups,
intermolecular sites, ionisation and an explicitly bounded partition model.
The middle demos use Newman/Fischer views, populations, ee, assumed rate laws,
E2 geometry, oxidation states, elementary mechanisms and stoichiometric limits.
The later demos cover amine protonation, resonance, spectral candidates,
TLC, short routes and conditional evidence claims.

## Implementation

Three parallel workers implemented the day groups. A shared interface supplies
labelled native controls, responsive SVG, status updates and numerical formatting.
The runtimes use classic local scripts, with no imports, fetch, CDN, build step,
student accounts, stored answers or external service. They work from the local
rendered files as well as the GitHub Pages directory structure.

The diagrams are teaching models, not molecular dynamics or operating protocols.
Supplied constants, simulated observations and selected pathways are identified
on their pages. In particular, extraction fixes aqueous pH and partitions only
neutral HA; rate models do not select a real substrate mechanism; qualitative
spectra do not prove purity; and the final reference-pattern example is not an
aspirin identification spectrum borrowed from another teaching sample.

## Review and verification

- 80 pure-model checks pass: 17 foundation, 29 reactivity and 34 application
  anchors/boundaries. They cover, among other things, hydrogen/formula counts,
  Fischer transformations, acid/base fractions, partition mass balance,
  population normalisation, charge conservation, limiting amounts, invalid
  masses, Rf bounds and evidence limitations.
- Workers independently checked the foundation and application chemistry.
  Parent review checked the intermediate models and their existing course
  assumptions. An exchangeable OH feature was broadened to match its label.
- All 25 rendered demos load in Chrome from local files with external HTTP(S)
  requests blocked. Each responds to real controls, has a functioning reset,
  and fits desktop/mobile widths. The checks include alternate selections,
  slider endpoints, blank numeric input and a keyboard-driven pH change.
- Every rendered sidebar contains the 25 direct demo destinations. Static
  checks verify page roots, local runtime paths, script syntax and built links.
- Visual inspection led to compact comparison layouts and removal of a clipped
  caption. Representative desktop/mobile ionisation and spectroscopy pages,
  as well as geometry views, were inspected after those changes.
- The complete site renders without unresolved-link warnings. Existing course,
  paper-bank, R example and clean-report checks still pass.

`make demo-check` runs the source/model checks. Optional `make demo-browser`
runs the local-browser checks and writes `demo-browser-check.json` here.
`make check` includes the demo source/model checks alongside the previous
course verification. Student-machine/classroom pilots remain useful for pace
and projection readability.
