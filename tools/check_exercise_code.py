"""Execute the supplied student R examples from a clean extracted project.

Restrict this to complete examples in exercise pages and the TLC reference;
setup installation commands and lecture fragments are deliberately not run.
"""
import re
import subprocess
import tempfile
from pathlib import Path
from zipfile import ZipFile

ROOT=Path(__file__).resolve().parents[1]
pages=sorted(ROOT.glob('*/exercise.qmd'))
count=0
with tempfile.TemporaryDirectory(prefix='orgchem-r-exercises-') as td:
    project=Path(td)
    with ZipFile(ROOT/'downloads/organic-chemistry-student.zip') as z: z.extractall(project)
    for page in pages:
        blocks=re.findall(r'^```r\n(.*?)^```',page.read_text(),re.S|re.M)
        if not blocks: continue
        # Disable interactive graphics; each example uses real bundled inputs.
        script=project/'exercise-check.R'
        script.write_text('pdf(NULL)\n'+ '\n\n'.join(blocks)+'\ndev.off()\n')
        result=subprocess.run(['Rscript','--vanilla',str(script)],cwd=project,
                              capture_output=True,text=True)
        if result.returncode:
            print(f'FAILED: {page.relative_to(ROOT)}\n{result.stdout}\n{result.stderr}')
            raise SystemExit(result.returncode)
        count+=1
        print(f'R example passed: {page.relative_to(ROOT)}')
print(f'All supplied complete R examples passed in {count} pages.')
