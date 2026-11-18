# Organic Chemistry 1 — proposed redesign

Working design, 18 November 2026. **Target: the new 5KE231, 7.5 credits**, for
the Life Science programme's 2027 intake. This is not an approved syllabus.
Evidence, textbook links and uncertainties are recorded in [sources](plan/sources.md).

**Implementation status, 19 November 2026:** the first teaching draft now
contains all 15 units, 45 unit pages plus seven course/reference pages,
150 original SVG diagrams, R plots, practice and explained self-checks,
constructed teaching data and a working student Rmd project. The pedagogical
proposal below remains the basis for review. Department-approved laboratory
protocols, final syllabus alignment and student pacing pilots remain open.

## Purpose and boundaries

Students should learn to read organic structures, reason about electron flow,
predict a small set of important transformations, and interpret experimental
evidence. R is their ordinary tool for calculations, tables and plots; R Markdown
is their ordinary tool for a reproducible laboratory report. Most chemistry
learning remains drawing, explaining, modelling and solving problems.

The local programme record lists Organic Chemistry II as a **year-3 elective**.
Consequently, this course must give every student a usable foundation for
biochemistry without depending on their taking II. Keep simple carbonyl and
acyl chemistry here, even if that requires reducing the breadth of hydrocarbon
reaction coverage. II can deepen synthesis, selectivity and reaction families.
That division is our proposal; II's detailed syllabus is unknown.

Use a small recurring set of molecules: ethanol/ethanal/ethanoic acid for
functional-group transformations; a simple ester for synthesis and hydrolysis;
lactic acid for stereochemistry; an amino acid for acid–base behaviour; aspirin
as an optional integrative case. Introduce biological relevance through the
chemical feature under discussion, supplying the necessary biological context.

## Proposed learning outcomes

By the end, students should be able to:

1. Translate between molecular, Lewis, condensed, skeletal and basic 3D
   representations; identify functional groups and assign formal charges.
2. Use polarity, resonance, intermolecular interactions and acid–base reasoning
   to explain selected properties and reactive sites.
3. Distinguish conformation from configuration, identify common stereochemical
   relationships and assign simple R/S and E/Z configurations.
4. Draw chemically valid elementary mechanisms, following electron pairs,
   charges, proton transfers and intermediates.
5. Predict and justify representative substitution, elimination, addition,
   carbonyl addition, acyl substitution and alcohol redox transformations.
6. Combine a small reaction toolkit into a justified one- or two-step route,
   including likely competing reactions and useful experimental checks.
7. Combine simple IR and one-dimensional NMR evidence with molecular formula
   and chemical reasoning to distinguish candidate structures.
8. Record and evaluate basic laboratory work, calculate stoichiometry and
   isolated yield, and communicate evidence and limitations in an Rmd report.

Outcomes 4–7 require substantial work on paper or a structure editor. Producing
an R plot is not evidence that a student can explain a mechanism.

## Content selection and depth

| Area | Required capability in Organic Chemistry 1 | Boundary for this draft |
| --- | --- | --- |
| Representations and naming | Read skeletal structures, recover implicit atoms, name simple examples and locate functional groups | Avoid an exhaustive nomenclature catalogue |
| Bonding and resonance | Use geometry, sigma/pi bonds, charge and valid resonance contributors to explain a reactive site | No orbital calculations or quantum-mechanical derivations |
| Acids, bases and properties | Compare stabilisation of conjugate bases; predict proton transfer; interpret supplied pKa values and aqueous ionisation | Teach model assumptions; no numerical equilibrium solver |
| Conformation | Newman projections, anti/gauche, basic cyclohexane chairs and axial/equatorial preference | No extensive substituted-ring energy calculations |
| Stereochemistry | Simple R/S, E/Z, enantiomers, diastereomers, meso examples and conformers | R/S does not predict optical rotation; biological effects need evidence |
| Mechanisms | Electron source → sink, nucleophile/electrophile, leaving group, intermediate versus transition state | Track atoms and charge in every worked mechanism |
| Substitution/elimination | SN2 inversion; SN1 intermediates; E2 geometry; recognise E1 and competing pathways | Use a few justified cases; avoid a memorised universal decision tree |
| Unsaturated molecules | Representative electrophilic additions, regioselectivity and an introductory stereochemical example | Alkynes at recognition level; reduce reagent variants |
| Alcohols and carbonyls | Relate alcohol oxidation/reduction to aldehydes/ketones; explain simple nucleophilic addition and a hemiacetal | Supply selected reagents; omit broad organometallic and protecting-group chemistry |
| Acids, esters and amides | Explain elementary acyl addition–elimination, ester formation/hydrolysis and amide stabilisation | Amide formation needs activation or appropriate conditions; not simply mixing acid and amine |
| Aromatics and amines | Recognise aromatic stabilisation and simple heterocycles; compare selected amine basicities | One illustrative aromatic substitution if time; directing-group catalogues deferred |
| Spectroscopy | IR functional groups; simple proton NMR shifts, integrations and splitting; basic carbon NMR environments | No 2D NMR or automated structure elucidation; mass spectra optional |
| Synthesis and practical work | One- or two-step transformations, separation/purity concepts, observations, yield and evidence | Use department-approved protocols; exact lab allocation remains open |

Candidates for II include fuller aromatic substitution and directing effects,
radical chemistry, organometallic reagents, enolates and condensations,
conjugate addition, advanced selectivity, longer synthesis and more demanding
spectral problems. A short recognition example may appear here without making
the whole family examinable. Detailed metabolism belongs in Biochemistry.

## Five-week teaching route

Assume five full-time weeks for planning. These are study blocks, not a confirmed
timetable. Each week has three short lecture decks, two tutorial/practical blocks
and protected independent work. Revisit representations, charge and electron flow
throughout. Introduce IR early, then reuse it when reactions change functional groups.

Each deck initially targets **16 slides**, including opening and recap slides,
with **at least 10 substantive visual slides**. Adjust slide count to teaching
pace while retaining the per-deck visual threshold described below.

| Week | Unit / proposed directory | Chemistry question and activity | Planned supporting visuals |
| --- | --- | --- | --- |
| 1 | 01 `structures/` | What does a drawn molecule tell us? Translate structures, count atoms, identify functional groups, use simple names | Representations of one molecule, annotated functional groups, geometry models |
| 1 | 02 `properties/` | Why do similar molecules behave differently? Polarity, hydrogen bonding, intermolecular forces and introductory IR | Interaction diagrams, real solvent/property examples, annotated IR comparisons |
| 1 | 03 `acids-bases/` | Where will a proton move? Resonance, conjugate bases, pKa and ionisation | Electron-flow diagrams, conjugate pairs, aqueous fraction curves |
| 2 | 04 `conformations/` | Is this a new molecule or a rotated view? Newman projections and chair conformations | Linked molecular views, rotation sequence, energy profile, chair flip |
| 2 | 05 `stereochemistry/` | Are these structures identical, mirror images or different stereoisomers? | Wedge/dash and model pairs, symmetry overlays, R/S and E/Z worked examples |
| 2 | 06 `mechanisms/` | How do we account for a reaction? Nucleophiles, electrophiles, electron flow and energy | Curved-arrow steps, charge tracking, intermediate/transition-state comparison |
| 3 | 07 `substitution/` | How can one group replace another? Contrast a few SN2/SN1 cases | Backside attack, inversion, carbocation geometry, qualitative rate comparisons |
| 3 | 08 `elimination-addition/` | When does a pi bond form or react? E2/E1 recognition and representative additions | Anti-periplanar geometry, substitution/elimination comparison, addition mechanisms |
| 3 | 09 `alcohols-redox/` | How are alcohols, aldehydes and ketones connected? | Transformation map, oxidation-level comparison, before/after IR |
| 4 | 10 `carbonyl-addition/` | Why is a carbonyl carbon attacked? Simple addition, proton transfers and hemiacetals | Polarised carbonyl, stepwise addition, tetrahedral product, ring/chain example |
| 4 | 11 `acyl-chemistry/` | Why do esters and amides behave differently? Elementary addition–elimination and hydrolysis | Tetrahedral intermediate, leaving-group comparison, ester/amide resonance |
| 4 | 12 `amines-aromatics/` | How do nitrogen and aromatic rings change properties? | Lone-pair availability, aromatic/nonaromatic contrasts, amino-acid charge forms |
| 5 | 13 `spectroscopy/` | Which candidate fits all the evidence? Integrate formula, IR and simple proton/carbon NMR | Authentic or labelled simulated spectra, candidate elimination, peak annotations |
| 5 | 14 `synthesis-evidence/` | Did we make the intended compound? One- or two-step synthesis, purity and yield | Reaction map, TLC plate, isolation/purification photos, spectra and yield table |
| 5 | 15 `integration/` | Can we explain an unfamiliar example using the same principles? | Annotated molecule, mechanism and evidence panels, student explanation prompts |

Suggested remaining blocks: week 1 representation practice and Rmd restart;
week 2 model-building and mixed stereochemistry/mechanism problems; week 3
mechanism clinic and first laboratory session; week 4 carbonyl problems and
second laboratory session; week 5 structure puzzles, report discussion and
revision. Pilot the week-4 load first. If overloaded, cut aromatic reaction
coverage and extra reagent variants before cutting core carbonyl reasoning.

## R as a routine chemistry tool

Reuse `readr`, `dplyr`, `ggplot2`, ordinary arithmetic and the familiar report
structure from `introcomputing`. No new package ecosystem is needed. Supply a
small project with `data/`, `output/` and `report.Rmd`, relative paths and an
untouched raw-data copy. Provide a 30–60 minute optional catch-up walkthrough
if the introductory course has not yet run or skills need refreshing.

Each task follows **predict → calculate/plot → check → explain the chemistry**.
Start with a manual example and reduce template guidance over time. Students
should leave with saved work they can reopen, not only console output.

| Placement | Task and chemical purpose | Familiar R work | Time / required output |
| --- | --- | --- | --- |
| Week 1 | Reopen a starter report and calculate amounts from mass and supplied molar masses | Named quantities, unit conversions, derived columns, Knit | 30 min; one checked calculation table |
| Week 1 | Plot the aqueous ionised fraction of a supplied monoprotic acid across pH | Supplied formula, `seq`, powers, labelled line plot | 40 min; curve plus explanations at pKa and pKa ± 1 |
| Week 2 | Compare two supplied conformer free energies at one temperature | `exp`, normalisation, bar plot | 35 min; populations and a prediction checked against the figure |
| Week 2 | Convert a simple two-enantiomer composition into enantiomeric excess | Table arithmetic and a composition plot | 25 min; explain why ee is not the percentage of the major enantiomer |
| Week 3 | Prepare a balanced-reaction amount table and identify the limiting reagent | CSV import, amount/stoichiometric-coefficient comparison | 45 min; theoretical product amount with units |
| Week 3–4 | Evaluate supplied or measured TLC distances and keep sample identity | Derived `Rf` column, simple table/plot | 30 min; chemical interpretation and a statement of what TLC cannot establish |
| Week 5 | Plot supplied IR traces and compare candidate functional groups | Read numeric columns, `geom_line`, reversed wavenumber axis | 45 min; annotated evidence supporting a candidate |
| Week 5 | Calculate isolated yield and assemble experimental evidence | Reuse amount table, measured mass, inline R results and figure captions | 50 min; draft results section and interpretation |

These eight tasks total five hours. Reserve approximately another three hours
for report assembly and troubleshooting: **about eight hours of explicit R/Rmd
work**, integrated into the chemistry workload. Chemistry problem solving gets
far more time. Optional kinetic plots must replace an activity if added.

Supply equations and assumptions explicitly. For example, the aqueous
monoprotic-acid model is `fraction_A_minus = 1 / (1 + 10^(pKa - pH))` using a
supplied pKa and the stated dilute-solution approximation. At pH = pKa the
fraction is 0.5; at pKa + 1 it is about 0.909. This is not a solvent partition
model. Do not infer extraction recovery from that curve alone.

For conformers supply relative **free energies**, temperature and consistent
units, using `exp(-delta_G / (R * temperature_K))` then normalising the weights;
state that any degeneracy is already accounted for. Keep this guided and omit
it if the mathematical load obscures conformation. For spectra, R presents the
observations; students assign structures and justify them themselves.

Do not introduce t-tests, machine learning, docking, quantum chemistry, spectral
classification, cheminformatics programming or database pipelines as core work.

## Practical work and reporting

Provisionally budget two six-hour laboratory sessions, subject to local
facilities and approved protocols. One could emphasise separation, TLC and
purification; another a simple synthesis followed by identification. Ester or
aspirin chemistry is a candidate context, not a selected or approved experiment.
Do not add reaction instructions until the department's protocols are available.

Students retain dated primary observations and sample identities in a lab
journal. The Rmd report draws on those records and saved data. It includes the
chemical question, reaction scheme, relevant conditions, amount table,
observations, calculations, spectra/TLC, interpretation and limitations.
An annotated photograph or exported structure drawing can accompany the code.

Distinguish isolated yield, conversion and purity. Repeated instrument readings
are not independent syntheses. A clean-looking TLC spot or plausible isolated
mass does not establish identity or purity. When student measurements are not
available, distribute explicitly labelled teaching data with units, provenance
and a data dictionary; synthetic spectra must be labelled as such.

For a checkable amount/yield exercise, supply a balanced 1:1 reaction with
10.0 mmol of limiting starting material, 15.0 mmol of the other reagent and
product molar mass 180.16 g/mol. Theoretical product mass is 1.8016 g; an isolated
1.20 g corresponds to about 66.6% yield, without establishing purity.

## Workload and checks of learning

Illustrative planning envelope: **200 hours**, to be reconciled with the local
credit-hour convention and actual timetable. Allocate 15 hours to lectures,
20 to tutorials, 8 to explicit R/Rmd work, 12 to laboratory attendance, 20 to
lab preparation/journalling, 85 to reading and independent chemistry problems,
and 40 to integration/revision/assessment. Hours are non-overlapping: do not
charge the R task again inside tutorial or report time.

Use short formative quizzes after each unit, combining structure recognition
with an explanation. Revisit prior concepts in later units. Examples include
repairing an invalid curved arrow, comparing conformers without changing
configuration, and rejecting a spectral candidate using two observations.

Proposed assessment evidence: individual structure/mechanism/spectroscopy
problems for outcomes 1–7; observed practical work, journal and report for
outcome 8; a short discussion that asks students to explain their own result.
Grade chemical correctness, evidence and reasoning, with reproducibility as
part of competent reporting. Do not make typing speed or decorative R code a
major determinant. Formal examination modes and weightings await the syllabus.

## Visual standard: more than half of every deck

The working target is **10 visual slides per 16-slide deck**, or about **150
of 240 slides (62.5%)** across 15 decks. The acceptance floor is **strictly
more than 50% in each deck**, counting all visible slides, including opening
and recap slides. Speaker notes are excluded. A whole-course average must not
conceal a text-heavy deck.

A qualifying slide has a meaningful structure drawing, mechanism diagram,
molecular model, spectrum, graph, experimental photo or other visual that
students must inspect to answer the chemical question. A decorative molecule,
logo or background photograph does not qualify. A small formula in a paragraph
does not qualify; an annotated reaction scheme or molecular representation can.

For an initial 16-slide deck, plan roughly six slides with structures,
mechanisms or spatial models, two with plots/spectra, and two with practical
images or other concept visuals. Adapt that mix to the topic; stereochemistry
will need more models, while spectroscopy will need more spectral evidence.
Use progressive states to teach electron flow, not to inflate the visual count.

Use editable chemical drawings for exact structures and mechanisms, verified
3D renders or model photographs for spatial relationships, and R for plots.
Raster illustrations can support nontechnical context, but exact chemistry
should be drawn and reviewed for valence, charge, arrow origin/destination,
stereochemistry and atom conservation. Do not rely on generated text or
generated molecular drawings for assessed chemical content.

Maintain an asset register with filename, slide use, source, licence/permission,
modifications, attribution and chemical review status. Prefer original drawings
and lab photos, or clearly reusable sources. Publisher images require suitable
reuse permission. Provide descriptive alt text, legible labels, high contrast
and cues beyond colour. Spectra need axes and units; photos need annotations
when the relevant observation is not obvious.

## Scaffold and implementation sequence

**1. Establish the minimal site.** Copy/adapt `_quarto.yml`, `styles.css`,
`theme.scss` and the topic structure from
`/home/mahogny/github/course/quantbio-3MB00A`. Keep Quarto HTML pages, Reveal.js
slides, search, docked sidebar, navigation controls, speaker notes and `docs/`
output. Rewrite course identity and navigation. The source Makefile's default
also stages files; replace that behaviour with a render-only default. Leave
unneeded demos, advanced course content and dependencies out of the scaffold.

**2. Pilot a complete unit.** Build `acids-bases/lec.qmd`, `quiz.qmd`,
`exercise.qmd`, its ionisation dataset/code and report section, plus at least
ten reviewed visual slides in the proposed 16-slide deck. This tests the blend
of chemistry, R, visuals and navigation. Check a second demanding example from
stereochemistry before scaling the visual approach.

**3. Write the chemistry spine.** Draft all 15 units with explicit outcomes,
selected textbook sections, worked problems, common errors and quizzes.
Approve chemical scope and prerequisites, particularly the week-4 load and
handover to Biochemistry. Keep optional material visibly separate.

**4. Build the exercises and practical package.** Add the eight short R tasks,
CSV dictionaries, a starter `.Rmd`, student ZIP, approved lab material when
available, and instructor solutions. Maintain continuity with the introductory
course's software and report conventions. Student reports stay `.Rmd` even
though the course website uses `.qmd`.

**5. Review and validate.** Review every chemical drawing and worked answer.
Render the site and inspect links, navigation, slide overflow and projected
readability. Count meaningful visual slides per deck, excluding speaker notes
and decorative images. Knit the starter report in a clean R session from an
extracted student download; check units and known numerical answers. Keep
instructor solutions, planning files and private assessment material outside
the public render and downloads. Pilot with beginning students and adjust time.

Suggested source layout, to create during implementation:

```text
_quarto.yml, Makefile, index.qmd, styles.css, theme.scss
structures/ ... integration/  # each: lec.qmd, quiz.qmd, exercise.qmd, figures
intro/                      # short R/Rmd restart and setup links
reference/                  # functional groups, vocabulary, reaction toolkit
data/                       # CSVs, dictionaries and provenance
templates/report.Rmd
downloads/                  # packaged student project
instructor/                 # solutions and teaching notes; excluded
plan/                       # evidence and asset/slide inventories; excluded
tools/                      # packaging and focused validation
docs/                       # rendered student site
```

The scaffold, pilot and full first teaching draft are implemented. The next
teaching milestone is syllabus reconciliation, local practical arrangements
and a student/projection pilot before classroom release.
