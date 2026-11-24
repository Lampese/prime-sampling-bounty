import json
from pathlib import Path
from functools import lru_cache
from fractions import Fraction as F
from math import comb,prod

from . import certificate as c

BASE=Path(__file__).resolve().parent
original=c.group
patterns=json.loads((BASE/'bank.json').read_text())
for r in range(6,22):
    patterns[str(r)]=[(h,n,t) for h in range(r+1) for n in range(h,r+1)
                      for t in (0,2,4) if t<=n-h]
@lru_cache(None)
def bank(R):
    delta, old, lower=original(R)
    opts=set(old)
    for h,n,t in patterns[str(R)]:
        assert 0<=h<=n<=R and t in (0,2,4,6,8)
        gh=prod(1-F(1,p) for p in c.ps[:h])
        qs=c.ps[h:n];dist=[F(1)]
        for q in qs:
            u=F(1,q)
            dist=[a*(1-u)+b*u for a,b in zip(dist+[0],[0]+dist)]
        vals=[1 if j==0 else (comb(j-1,t) if j>t else 0) for j in range(len(dist))]
        mu=gh*sum(w*f for w,f in zip(dist,vals))/delta
        var=gh*sum(w*f*f for w,f in zip(dist,vals))/delta**2-mu**2
        norm=2**h*gh/delta*sum(comb(n-h,i) for i in range(min(n-h,t)+1))
        if sum(F(1,p) for p in qs)<=1:
            e=[F(1)]+[F(0)]*t
            for q in qs:
                for j in range(t,0,-1):
                    e[j]+=(1-F(2,q))*e[j-1]
            norm=min(norm,2**h*gh*sum(e)/delta-mu)
        opts.add((c.up(mu-1),c.up(var),c.up(norm),h,t))
    kept=[]
    for b in sorted(opts):
        if not any(a[0]<=b[0] and a[1]<=b[1] and a[2]<=b[2] for a in kept):kept.append(b)
    return delta,kept,lower
