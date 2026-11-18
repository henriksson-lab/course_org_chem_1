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

The draft includes 15 illustrated lecture decks, worked practice, explained
self-checks, eight short R tasks, constructed datasets and dictionaries, a
working Rmd report and a packaged student project. Outcomes, teaching hours
and assessment suggestions remain design choices, not an approved syllabus.
Actual laboratory protocols and local arrangements must be supplied by the
department. The provided practical work covers preparation and interpretation.

Run `make render` (or `make`) to package the download and render all 52 student
pages into `docs/`; `make preview` opens a local preview. `make check` verifies
the course structure, figures, links, downloads, complete R exercise examples
and knitting the extracted student report. `make figures` regenerates all
original figures and their register. R, Quarto and Python are required; see
the instructor guide for package details. No build target stages or commits.

Planning documents and instructor tools are excluded from the rendered site
and student ZIP. The source and built course remain in this repository.

Design requirements: a self-contained foundation for students who will not
choose Organic Chemistry II; chemistry reasoning before calculation;
student reports in `.Rmd`; meaningful supporting visuals on more than half
of the slides in every deck.
