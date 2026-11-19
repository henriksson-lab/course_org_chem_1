# Brief: rework one lecture deck (instructor planning note)

The lecture decks (`<unit>/lec.qmd`, Quarto revealjs) were generated to fill a
40-seconds-per-slide quota. They have about 55 slides each; one-line
"Predict:"/"Resolve:" slide pairs; `figures/NN-pace.svg` "figures" that are
only text in a box; copy-pasted boilerplate speaker notes; no opening or
closing slide; and a rubric-like tone. Nobody wants to listen to them.

The model to follow is **`structures/lec.qmd`**, which has already been
reworked. Read it fully before starting, together with its figure script
`tools/make_structure_figures.py`. The house style comes from the quantbio
course, e.g. `/data/henriksson/github/course/quantbio-3MB00A/pca/lec.qmd`.

## Required deck shape

- **Front matter:** copy it from `structures/lec.qmd` (`author: Organic Chemistry 1`,
  `history: false`, theme and size settings). Drop `session-minutes` and
  `pacing-seconds-per-slide`. Set the subtitle to `Week N · Day D — <short hook>`,
  keeping the week from the old deck.
- **`## About this lecture`:** open with a concrete, interesting motivating
  case (biology, medicine, the lab, everyday life). Then give 4–5 short
  bullets on what we learn, and one line linking back to the previous lecture.
- **15–19 content slides** in a clear narrative order. Each slide covers one
  idea and has a real figure or a table, plus 1–3 sentences or 2–4 short
  bullets. Use worked examples. Avoid near-empty slides.
- **"Try it" boxes:** fold questions into the slide they belong to, using
  exactly this markup:

  ```
  ::: {.try}
  **Try it:** question

  ::: {.fragment}
  **Answer.** One-sentence explanation.
  :::
  :::
  ```

  Use about 6–9 per deck. Never put a question and its answer on separate slides.
- **`## Mistakes we see every year`** (optional but encouraged): a table with
  the columns "mistake | fix".
- **`## Summary`:** a "concept | key fact" table of at most about 9 short rows,
  then a single line: `**Next:** … (Day X). **Practice:** [exercises](exercise.qmd) · [self-check](quiz.qmd)`.
  For the last unit, point to revision instead of a next lecture.
- **Speaker notes** (`::: notes`) on every slide: specific, conversational
  guidance on what to say and which misconceptions to probe. Never use generic
  boilerplate.
- **Style:** plain, friendly English with contractions. No em-dash runs and
  no rubric/audit jargon ("audit", "bring these tools into the room"). Mention
  connections to other days ("comes back in Week 4").
- **Content coverage:** keep everything the old deck taught that the unit's
  `exercise.qmd` and `quiz.qmd` rely on. Read those to know what students
  must be able to do. The chemistry must be correct, so double-check every
  formula, H count, charge, pKa, name and stereodescriptor.

## Figures

- Draw molecules as real skeletal structures with RDKit through
  `tools/molfig.py` (read it). It provides:
  - `place(fig, smiles, x, y, numbers=, highlight_atoms=, explicit_h=, bond_px=, font_px=, rotate=, flip=)`,
    which returns the atom coordinates;
  - `curved_arrow(fig, start, end, bend=)` for electron-pushing arrows;
  - `lone_pair(fig, atom_xy, angle)`.
- Canvas and text come from `tools/chemfig.py` (`Figure`, `text`, `line`, `rect`, …).
  Stereo is drawn from SMILES `@`/`@@` (wedges). Newman/Fischer projections,
  energy diagrams, plots and spectra are not RDKit's job: reuse the existing
  good figures in `<unit>/figures/` (those not named `*-pace.svg`), or draw
  them with chemfig.
- Write your figure script as **`tools/make_sk_<unit>.py`**. Model it on
  `tools/make_structure_figures.py`, write to `<unit>/figures/sk-*.svg`, and
  run it.
- Do **not** edit shared files: `tools/molfig.py`, `tools/chemfig.py`,
  `theme.scss`, `Makefile`, `_quarto.yml`, `tools/check_course.py`, or any
  other unit. If molfig lacks something, work around it in your own script
  and mention it in your report.
- Use 960-px-wide canvases, typically 330–470 px tall. Molecules need to be
  big enough to read on a projector (bond_px about 50–80).
- No `*-pace.svg` figures and no "figures" that are just text in a box.
  Every image needs descriptive non-empty alt text: `![Title](figures/x.svg){width=..% fig-alt="..."}`.

## Verify (mandatory)

1. Render: `XDG_CACHE_HOME=$PWD/.quarto/cache quarto render <unit>/lec.qmd`
2. Screenshot and check for overflow:
   `python3 tools/slide_screenshots.py docs/<unit>/lec.html <scratch-dir>`.
   It prints `OVERFLOW` for any slide whose content runs off the slide. Fix
   every overflow (shorter text, smaller figure width, or a split slide).
3. Look at the screenshots yourself (montage them with ImageMagick
   `montage a.png b.png -tile 3x3 -geometry 640x360+4+4 m.png` and Read the
   image). Check that text is readable, figures are legible and nothing
   overlaps. Iterate until the deck looks good.
4. Run `python3 tools/check_course.py`; it must end with "pass".

## Deliver

Commit only `<unit>/lec.qmd`, `<unit>/figures/sk-*.svg` and
`tools/make_sk_<unit>.py` in your worktree. Do not commit `docs/` or
`plan/*.csv`. Leave the old `*-pace.svg` files in place; cleanup happens
centrally. Report in under 150 words: slide count, figures made, anything
you were unsure about chemically, and any molfig limitations you hit.
