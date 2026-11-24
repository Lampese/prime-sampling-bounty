from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

from .sieved_det import install as determinant

def install(f):
    from ..spectrum import squarerates as sr
    from . import allsquarecap as av

    determinant(f)
    ceil=av.ceildiv
    @lru_cache(maxsize=20000)
    def cap(a,b,v,p0,p1):
        if b<=a or v<=0:return v+1
        assert a>194 and v<=b
        T=max(1,isqrt(p1//b));assert T<p0
        G1=max(b*T,ceil(p1,T))
        G0=max(b,min(a*a,ceil(p0,b+v)))
        if G0>G1:return 0
        left=av.row_bound(a,b,v,p0,G0)
        right=av.row_bound(a,b,v,p0,G1)
        best=min(v+1,min(max(x,y) for x,y in zip(left,right)))
        slope=F((2*b+v)*v,2*b*p0);jmax=slope*G1
        for shift in (2,4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,80):
            g=max(G0,G1>>shift)
            low=min(max(x,y) for x,y in zip(left,av.row_bound(a,b,v,p0,g))) if g>G0 else 0
            if low>=best:continue
            kmax=F(v*b,g)+jmax/2
            area=slope*v*b+jmax*jmax/2+3*(jmax+kmax)+9
            high=av.bill((jmax+3).__ceil__(),(kmax+3).__ceil__(),area.__ceil__(),b,p0)
            if high is not None:best=min(best,max(low,high))
        return best
    prev=sr.square_cap
    def joined(*args):return min(prev(*args),cap(*args))
    sr.square_cap=joined
