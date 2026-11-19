# Instructor notes and review guide

This is the first teaching draft for the new Organic Chemistry 1, not a formal
syllabus. Planning and instructor files are excluded from the student website
by an explicit list of rendered student pages. The student ZIP uses an allowlist.

## Before teaching

Reconcile the draft with the new 5KE231 syllabus, actual reading list, timetable
and examination rules. Confirm the relationship to Fundamentals of Chemistry,
Data Analysis, Biochemistry and the elective Organic Chemistry II. Supply the
department-approved practical protocols and local arrangements. The current
site supports preparation/interpretation and does not prescribe lab operations.

Pilot student pacing, especially stereochemistry and the week-4 carbonyl load.
Use physical molecular models where available. Keep spatial drawing and
electron-flow practice central; reduce reagent variants before reducing that
practice. Guide the optional conformer free-energy model if its mathematics
distracts from the chemistry.

## Session pacing and daily paper practice

Each of the 25 study days also has a short interactive exploration, linked
directly from its day in the left sidebar to the Exercises `#explore` section. Use its 5–10 minute prediction/change/explanation
task within an existing worked-example or drawing block, or as optional
independent exploration. It does not add another compulsory block to each
45-minute lecture. The daily paper problems and R report remain important.

All workbenches are local classic JavaScript and SVG, with no installation,
account, CDN or data submission. `demos/index.qmd` lists their purposes;
`plan/demos.json` maps exercise anchors, shared partials, compatibility URLs, runtimes and days. `make demo-check` checks the
local runtimes, sidebar destinations and pure chemical/numerical anchors.
Optional `make demo-browser` opens every demo in Chrome with external requests
blocked, changes real controls, checks reset and invalid inputs, and inspects
wide/narrow layouts. Its results are in `plan/demo-browser-check.json`.

Each lecture targets 45 minutes, including its brief prediction, drawing and
discussion tasks. Around 40 seconds per slide is a useful planning estimate,
not a slide quota. Retain a slide only if it contributes an explanation, a
worked operation, a comparison or a meaningful check. Natural slide counts
vary by topic. Use speaker notes to plan the longer worked steps and pauses.

The 25-day route is in the sidebar (`_quarto.yml`); review days are `practice/dayNN.qmd`. Each lecture unit has at least
20 non-R problems, and ten consolidation days have at least 12 each: a bank
of at least 420 paper problems, separate from quizzes and R activities.
Core/transfer/challenge labels keep daily required time realistic while
providing substantial independent practice and revision. Every problem should
have all needed inputs, a targeted hint and a reasoned worked solution.
Students should be able to use this route without the textbook.

Browser Back should leave a presentation in one action. Shared `history: false`
restores the source scaffold's behaviour; slide navigation may update the
current URL but must not push additional browser history entries. The optional
browser check reproduces entry from the homepage, several slide changes and
one Back action.

## Solutions and teaching support

Each unit's practice and self-check pages contain collapsible worked checks.
These are intentionally student-visible formative solutions. Lecture speaker
notes contain additional teaching prompts and common errors. There are no
private exam questions or promised examination weightings in this repository.

The fully working report is a scaffold, not a completed scientific conclusion.
Its interpretation prompts remain for students. `check-calculations.R` records
independent numerical anchors; the expected limiting amounts are 10 and 15 mmol,
the theoretical product mass is 1.8016 g and the simulated isolated yield is
66.6%. TLC values are 0.24 and 0.64, and the crude lane has two distinct spots.

## Data and scientific limits

All data are constructed. The IR traces are schematic Gaussian bands; NMR is
a simplified first-order peak list. Neither is a substitute for real spectral
examples in a taught laboratory. The starter's ethanol/ethyl-ethanoate IR
comparison is unrelated to the aspirin isolation dataset. Preserve sample
identity in student discussions.

The conformer table supplies a two-state free-energy model with the combined
gauche state's degeneracy included. Do not add a second degeneracy factor.
Approximate aqueous pKa values ignore differences in measurement conditions;
the ionisation calculation is not a solvent extraction or titration model.

## Rebuild and validation

`make render` packages the student download and renders the website. `make`
does the same and never stages or commits files. `make figures` regenerates
the pilot and R plots; the other figure generators are listed in the Makefile.
`python3 tools/build_data.py` regenerates the fixed teaching tables.

After rendering, `make check` checks source completeness, SVG XML and alt text,
per-deck image quotas, built local links, exclusion boundaries, numerical
answers and knitting from an extracted student ZIP in a new R process.
The slide inventory is saved in `plan/slide-inventory.csv`.

Optional `make layout` uses installed Chrome and Python's `websocket-client`
to open only the local rendered decks, load every figure, and check all slide
bounds and browser Back behaviour. It writes `plan/layout-check.json`,
`plan/navigation-check.json` and representative screenshots under
`/tmp/`. This check needs the local browser to be runnable; in a restricted
executor it may need native execution approval. It is separate from `make check`.

The image quota check is a proxy, not a chemical or visual review. Inspect
valence, formal charges, atom conservation, arrow origins/destinations,
stereochemical views, spectral axes and legibility. Check representative
slides on the actual teaching projector. Text-only cards do not satisfy the
meaningful-visual requirement even if saved as an image.

The Makefile puts Quarto's compilation cache inside the ignored
`.quarto/cache/` directory so builds also work in a restricted workspace.
Override `QUARTO_CACHE_ROOT` if needed. Direct `quarto render` invocations use
Quarto's normal cache location unless you configure `XDG_CACHE_HOME` yourself.
