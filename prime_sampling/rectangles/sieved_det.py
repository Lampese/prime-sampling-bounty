from math import comb,isqrt
def balanced(n,f):
 if not f:return 10**9 if n else 0
 d,r=divmod(n,f);return (f-r)*comb(d,2)+r*comb(d+1,2)
def gain(k):
 m=2*k+1
 ret=1
 for p in range(2,2*m+1):
  if any(not p%d for d in range(2,isqrt(p)+1)):continue
  f=p-1
  while f<m:
   cost=min([balanced(m,f)]+[(s-k)*(s-k-1)//2+balanced(m-s,f-1) for s in range(m+1)])
   assert cost>=0
   ret*=p**cost;f*=p
 return ret
from functools import lru_cache
@lru_cache(None)
def alpha(n):
    assert n in (7,13,25,49,97)
    k=(n-1)//2;power=n*(n-1);weight=n**n*2**(2*(n-1))
    budget=400**power*gain(k)**2
    lo,hi=100,400
    while hi**power*weight<=budget:lo,hi=hi,2*hi
    while hi-lo>1:
        m=(hi+lo)//2
        if m**power*weight<=budget:lo=m
        else:hi=m
    return lo
def install(f):
    from . import allsquarecap, lattice

    lattice.det_factor=alpha
    allsquarecap.det_factor=alpha
    old=f.rates
    def rates(a,*args):
        assert a>194
        return old(a,*args)
    f.rates=rates
