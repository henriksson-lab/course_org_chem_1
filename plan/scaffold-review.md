# Quantbio scaffold review

Compared against `/home/mahogny/github/course/quantbio-3MB00A`:
`_quarto.yml`, `designprinciples.md`, `plan/course-order.md`, `Makefile`,
`styles.css`, `solving-systems/lec.qmd`, `solving-systems/exercise.qmd`
and `demos/README.md`.

## Conventions retained

- Quarto website rendered into tracked `docs/`, with a docked, searchable sidebar.
  Sidebar order defines the teaching route; stable topic folders hold `lec.qmd`,
  `exercise.qmd` and `quiz.qmd`.
- Lecture, Exercises, self-check and orange Demo links sit together. Demos are
  embedded in Exercises, with sidebar links directly to `#explore`.
- Collapsible hints and worked explanations support an attempt before revealing
  an answer. Chemistry drawings and small cases precede computation.
- Reveal slides use 1280 × 720, visible controls, menu, no transitions and
  `history: false`. Browser Back returns to the course home.
- Demos use local classic JavaScript and SVG without runtime network dependencies,
  including when opening the rendered course from local files.
- Lecture links are green, quiz links blue, demo anchors orange and Exercises
  italic. Week colours preserve the scaffold's visual grouping.

## Deliberate adaptations

- Five weeks contain five Day sections each, as requested. Review days reuse
  relevant lecture links rather than adding decks for navigation alone.
- Chemistry lectures support 45-minute sessions. Quantbio's shorter recorded-video
  decks are a different teaching format; slide counts here follow the chemistry
  explanations and activities, with meaningful visuals on over half the slides.
- R and Rmd support selected chemistry tasks and laboratory reporting. Quantbio's
  rule pairing every mathematical formula with R does not apply to every organic
  chemistry operation; drawing, mechanisms and evidence remain central.
- Plain JavaScript needs no TypeScript bundling step. Shared demo partials live in
  `demos/_dayNN.qmd` and are included by each day's Exercises. Former public demo
  pages redirect to the exercise anchor, preserving existing links.
- `make render` builds without staging or committing. Publication remains an
  explicit Git action rather than a side effect of a preview/build.

## Verification

`make check` checks the 25 day groups, exercise/demo anchors, local targets,
chemistry problem counts, assets, offline scripts, pure models and R examples.
`make demo-browser` exercises every embedded demo offline at desktop and mobile
widths and checks representative old-URL redirects. `make layout` checks slide
bounds and browser Back behavior.

All three checks passed after rendering this restructure: 420 paper problems,
80 pure model checks, complete R examples and a clean-session student report;
25 embedded demos at both widths, including redirects for days 1, 15 and 25;
all 784 slides with no overflow or missing figures; representative lecture Back
checks returning to the home page. Desktop and mobile screenshots were also
visually reviewed for the day-grouped sidebar and embedded demo layout.
