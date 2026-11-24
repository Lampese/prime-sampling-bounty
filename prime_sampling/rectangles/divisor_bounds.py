from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt
from collections import Counter
PRIMES=[p for p in range(2,200) if all(p%d for d in range(2,isqrt(p)+1))]
@lru_cache(None)
def image(m):
    counts=Counter({1:1}); rem=m
    for p in PRIMES:
        if rem%p:continue
        rem//=p; new=Counter()
        for mult,n in counts.items():
            new[mult]+=n*(p-1); new[2*mult]+=n*((p-1)*(p-2)//2)
        counts=new
    assert rem==1
    return counts

@lru_cache(None)
def divisor_caps():
    small=[p for p in range(2,512) if all(p%d for d in range(2,isqrt(p)+1))]
    coeff={}
    for d in range(4,10):
        c=Q(1)
        for p in small:
            if p>=2**d:break
            j=0
            while (j+2)**d>p*(j+1)**d:j+=1
            c*=Q((j+1)**d,p**j)
        coeff[d]=c
    rows=[]
    for b in range(96,385,8):
        cap=2**b
        for d,c in coeff.items():
            n=(c*2**b).__ceil__(); lo,hi=0,1<<((n.bit_length()+d-1)//d)
            while hi-lo>1:
                v=(lo+hi)//2
                if v**d>=n:hi=v
                else:lo=v
            cap=min(cap,2*hi)
        rows.append((cap,[(2*(2*k+1),1<<(b*k//(4*k+2))) for k in (6,16,32)]))
    return rows
