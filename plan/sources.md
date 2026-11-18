# Evidence for the Organic Chemistry 1 redesign

Checked 18 November 2026. Distinguish published facts, earlier local research,
and our proposed teaching design.

## Programme and course evidence

| Source | What it supports | Limitation |
| --- | --- | --- |
| [Local programme overview](../../molbiol_program/README.md) and [programme record](../../molbiol_program/PROGRAMME.md) | The 2027 intake has compulsory Organic Chemistry 1, 5KE231, 7.5 credits, in year 1; Organic Chemistry II is a recommended chemistry elective, 7.5 credits, in year 3 | These are locally researched records, not official syllabi |
| [Official 2027 programme syllabus](https://www.umu.se/en/education/syllabus/nglsc-id4272172/) | Original source recorded by the programme research | Could not retrieve this page during this session; retain that qualification when citing the local findings |
| [Local Organic Chemistry 1 dossier](../../molbiol_program/courses/5KE231-organic-chemistry-1.md) | Target-course evidence and explicit uncertainty | No retrieved 5KE231 syllabus, reading list, formal outcomes, examination rules or lab allocation |
| [Official Organic Chemistry 5KE180 syllabus](https://www.umu.se/en/education/syllabus/5ke180/) | A 15-credit predecessor centred on structure, stereochemistry, mechanisms, synthesis and IR/NMR, with experimental work and reporting | Retrieved during this session; its workload and assessment do not define the new course |
| [Official Organic Chemistry for Engineers 5KE202 syllabus and literature](https://www.umu.se/utbildning/kurs-och-utbildningsplan/5ke202/) | Mandatory Solomons, Fryhle and Snyder, *Solomons' Organic Chemistry*, 12th Global edition, 2016, ISBN 9781119248972; reading list valid from 1 September 2025 | Retrieved during this session; adoption in this separate course does not confirm adoption in 5KE231 |
| [Local elective dossier](../../molbiol_program/courses/2027-chemistry-electives-organic-medicinal-inorganic.md) | Organic Chemistry II's programme position and unresolved 5KExxx code | Its proposed content is an estimate, not a retrieved course syllabus |

Exact-code searches did not establish a public 5KE231 syllabus. Direct English
and Swedish syllabus URLs could not be retrieved. This is an evidence gap,
not proof that no syllabus exists. Likewise, listing Organic Chemistry II in
the programme does not establish a particular offering date or its prerequisites.

## Textbook depth

Use **Solomons 12th Global edition as a provisional reference**, based on local
adoption evidence. It is not yet a student purchasing instruction.

The [Wiley Global edition contents](https://www.wiley.com/en-ie/shop/general-chemistry/solomons%27-organic-chemistry-global-edition-12th-edition-p-9781119248927)
were available in search results; direct retrieval returned 403. Wiley's
[12th-edition companion contents](https://bcs.wiley.com/he-bcs/Books?action=contents&bcsId=10134&itemId=1118875761)
also establish the broad progression. Our selections below are a proposed
route through the book, not evidence of Umeå's assigned chapters. Confirm
section numbering against the actual Global edition before publishing readings.
Regional editions can differ, particularly in later chapters.

Proposed reference map:

- Chapters 1–3: molecular representations, functional groups, IR and acid–base reasoning.
- Chapters 4–5: conformations and stereochemistry.
- Chapters 6–8: substitution, elimination and addition.
- Chapters 9, 11–12: introductory structure determination, alcohols and redox.
- Selected chapters 14, 16–17 and 20: aromaticity, carbonyl and acyl chemistry, amines.
- Selected examples from chapters 22–24: biological molecules as applications.

Do not assign all of those chapters in full. In particular, a life-science
foundation needs selected carbonyl chemistry even though the book places it
well after the first reaction chapters. This motivates reordering rather than
stopping at the middle of the textbook. Advanced enolate chemistry and a broad
synthetic reagent inventory are candidates for a later elective.

## Local teaching materials inspected

- [`introcomputing/contentplan.md`](../../introcomputing/contentplan.md): files,
  tables, R calculations, `readr`, `dplyr`, `ggplot2`, and reproducible reports.
- [`introcomputing/report.Rmd`](../../introcomputing/report.Rmd): familiar
  import → derived columns → table/plot → interpretation → Knit workflow.
- [`introcomputing/instructor/chemistry-context.md`](../../introcomputing/instructor/chemistry-context.md):
  chemistry examples proposed there are not necessarily implemented lessons.
- `/home/mahogny/github/course/quantbio-3MB00A`: inspected `_quarto.yml`,
  `Makefile`, `README.md`, `index.qmd`, `theme.scss`, `styles.css` and a lecture
  source. Suitable conventions include one topic directory, Reveal.js slides,
  exercises, quizzes, speaker notes, docked navigation and `docs/` output.

The programme research does not confirm that Data Analysis precedes Organic
Chemistry 1. Familiar R skills are therefore a desired baseline, with a supplied
starter report and a short catch-up route until sequencing is settled.

## Evidence to reconcile before final teaching release

Obtain the new course syllabus and reading list when available, coordinate
minimum entry chemistry with Fundamentals of Chemistry, and agree the handover
to Biochemistry and Organic Chemistry II. Confirm actual laboratory facilities,
approved experiments, teaching language, examination requirements and course
sequence. None of these gaps prevents drafting the proposed core.
