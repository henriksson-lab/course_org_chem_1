# Organic Chemistry 1 — teaching materials

First teaching draft for the **new 7.5-credit course, 5KE231**, in Umeå's
Life Science programme for the autumn 2027 intake.

- [Proposed content and implementation plan](contentplan.md)
- [Research evidence and remaining gaps](plan/sources.md)
- [Rendered student course](docs/index.html)
- [Instructor notes and checks](instructor/README.md)

The course keeps organic chemistry central and uses short R exercises and
R Markdown laboratory reports to sustain skills introduced in
[`introcomputing`](../introcomputing/README.md). The website uses the Quarto
scaffold from `/home/mahogny/github/course/quantbio-3MB00A`.

The draft includes 15 illustrated 45-minute lecture sessions, a 25-day practice
route with at least 420 non-R chemistry problems, explained self-checks, eight
short R tasks, 25 interactive chemistry demos embedded in daily Exercises with direct sidebar links,
constructed datasets and dictionaries, a working Rmd report and
a packaged student project. Every paper problem supplies its inputs, hint and
worked explanation; the local course can be used without the textbook.
Slide counts follow teaching needs; approximately 40 seconds per slide is
pacing guidance rather than a target. Outcomes, teaching hours
and assessment suggestions remain design choices, not an approved syllabus.
Actual laboratory protocols and local arrangements must be supplied by the
department. The provided practical work covers preparation and interpretation.

Run `make render` (or `make`) to package the download and render the student
site into `docs/`; `make preview` opens a local preview. `make check` verifies
the course structure, figures, links, downloads, complete R exercise examples
and knitting the extracted student report. `make figures` regenerates all
original figures and their register. Optional `make layout` checks every slide
in Chrome and confirms that browser Back returns from a lecture to the course
home. `make demo-check` checks all embedded daily demos and pure chemical models;
optional `make demo-browser` exercises controls and resets offline in Chrome
at desktop/mobile widths. The demos need no build step or external libraries.
R, Quarto, Python and Node.js (for checks) are required; see
the instructor guide for package details. No build target stages or commits.

Planning documents and instructor tools are excluded from the rendered site
and student ZIP. The source and built course remain in this repository.

Design requirements: a self-contained foundation for students who will not
choose Organic Chemistry II; chemistry reasoning before calculation;
student reports in `.Rmd`; meaningful supporting visuals on more than half
of the slides in every deck.

The sidebar follows Week → Day → Lecture, Exercises, Self-check and Demo.
Review days link to existing lectures. Demo links open the exercise’s `#explore`
section, following the quantbio scaffold; old demo URLs redirect there.
See [scaffold comparison](plan/scaffold-review.md) for the source conventions
and the deliberate chemistry adaptations.
