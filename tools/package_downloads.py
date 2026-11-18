"""Package only explicit student inputs; no planning or instructor files."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

ROOT=Path(__file__).resolve().parents[1]
dest=ROOT/'downloads'/'organic-chemistry-student.zip'
dest.parent.mkdir(exist_ok=True)
readme='''Organic Chemistry 1 — student project

1. Extract this ZIP into an ordinary folder.
2. Open organic-chemistry.Rproj in RStudio.
3. Open report.Rmd, replace the author, save and Knit to HTML.
4. Read data/README.md for units, provenance and sample meanings.

All data are constructed teaching examples. No experiment has been performed.
The IR traces are schematic models, not measured reference spectra. The IR
comparison in the report is unrelated to the aspirin isolation table.

Packages: readr, dplyr, ggplot2, knitr, rmarkdown.
Keep report.Rmd beside data/; keep raw CSVs unchanged.
Add your unit exercises to the report and save derived files in output/.
Restart R and Knit before sharing the report with its data and source.
'''
with ZipFile(dest,'w',ZIP_DEFLATED) as z:
    for f in sorted((ROOT/'data').iterdir()):
        if f.is_file() and f.suffix in ['.csv','.md']:
            z.write(f,'data/'+f.name)
    z.write(ROOT/'templates/report.Rmd','report.Rmd')
    z.write(ROOT/'organic-chemistry.Rproj','organic-chemistry.Rproj')
    z.writestr('README.txt',readme)
    z.writestr('output/.gitkeep','')
print(f'Packaged {dest.relative_to(ROOT)}')
