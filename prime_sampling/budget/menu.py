from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import json

from .certificate import up
from .tight import bank

here=Path(__file__).resolve().parent
source=json.loads((here/'rows.json').read_text())
@lru_cache(None)
def group(r):
    delta,opts,low=bank(r);opts=list(opts)
    for t,m,n,q,v in source:
        if t<=r:
            scale=10**10;opts.append((up(F(m,scale)/delta-1),up(F(v,scale)/delta**2),up(F(n,scale)/delta),0,0))
    kept=[]
    for row in sorted(set(opts)):
        if not any(all(x<=y for x,y in zip(k[:3],row[:3])) for k in kept):kept.append(row)
    return delta,kept,low
