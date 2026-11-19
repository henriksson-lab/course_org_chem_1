"""Expand the local Quarto includes used to share exercise/demo content."""
import re
from pathlib import Path

def expand_includes(path, stack=()):
    path=Path(path).resolve()
    if path in stack:raise ValueError('Circular source include: '+str(path))
    text=path.read_text()
    def include(match):
        target=match.group(1).strip('"\'')
        return expand_includes(path.parent/target,stack+(path,))
    return re.sub(r'\{\{<\s*include\s+(\S+)\s*>\}\}',include,text)
