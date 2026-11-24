from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

@lru_cache(None)
def wheel_expansion(W):
    terms = [(1,1)]
    den = 1
    n=W
    q=2
    while q*q <= n or n>1:
        if q*q>n:
            q=n
        if n%q==0:
            n//=q
            assert n%q
            den*=q-1
            terms=[(d*w,k*c) for d,k in terms for w,c in ((1,q-2),(q,1)) if k*c]
        q+=1
    assert sum(k for d,k in terms)==den
    return terms,den

def weighted_near(a,b,V,pmin,pmax,W,H):
    assert 2<=a<=b and 1<=V<=b and 2*W*H<pmin
    K=W*H
    numerator=0
    terms,den=wheel_expansion(W)
    for d, mult in terms:
        if K<d:
            continue
        candidates=[(V+1)*3*H]
        for aux in ((1,2) if a>2 else (1,)):
            D=isqrt((V*pmin*d)//(2*b*K*aux))
            if D<1:
                continue
            U=b*((pmax+D-1)//D + ((D+1)*aux*K+d-1)//d)
            f=0;power=a
            while power<=U:
                f+=1;power*=a
            assert f>0
            root=isqrt((2*V*b*K)//(pmin*d*aux))+1
            candidates.append(3*(H+1)*(f-1)+f*(4*H+6)*root)
        numerator+=mult*min(candidates)
    return F(numerator,den*3*H)
