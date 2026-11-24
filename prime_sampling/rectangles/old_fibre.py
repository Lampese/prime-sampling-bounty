from functools import lru_cache
from math import isqrt

from . import divisor_bounds as cent
def arc(width):
    return max(2,min(max(cap,min(a*(1+(width+d-1)//d) for a,d in samples))
                     for cap,samples in cent.divisor_caps()))
arc=lru_cache(16384)(arc)


def ce(a,b):return -(-a//b)
@lru_cache(None)
def wheels():
    m=1; vals=[]
    for p in cent.PRIMES[:14]:
        m*=p; image=cent.image(m)
        vals.append((m,tuple(image.items()),sum(image.values())))
    return vals
@lru_cache(maxsize=16384)
def cap(v,b,pm):
    A=8*b*v*v;line=arc(v);best=max(line,2+ce(A,pm))
    for m,im,count in wheels():
        den=pm*m*m; num=2*den+A
        if m<=510510:
            per=arc(ce(v,m))
            best=min(best,ce(sum(n*max(num,den*min(line,mult*per)) for mult,n in im),den))
        s=1+isqrt((A*line)//(pm*m*count))
        tn,td=(num,den) if num>s*den else (s,1)
        sub=pm*m*(tn-td)
        best=min(best,max(ce(count*tn,td),ce(line*(sub+A*td),sub)))
        if 2*count>=best:break
    return best
