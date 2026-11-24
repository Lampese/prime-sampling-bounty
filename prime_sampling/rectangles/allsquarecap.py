from functools import lru_cache
from math import isqrt
from fractions import Fraction as F

from .lattice import ceildiv,grid_power,det_factor
from .fine_square import optimized
from .allmods import upper as all_moduli

@lru_cache(None)
def prime_count(j):
    if j<20:return 0
    N=1<<j
    lo=F(693147180559,10**12);hi=F(693147180561,10**12)
    return (2*N/((j+1)*hi)-N/(j*lo-F('1.1'))).__ceil__()

def pow_lower(g,k):
    j=g.bit_length()-1
    return (grid_power(j,k)*((2<<j)-g)+grid_power(j+1,k)*(g-(1<<j)))//(1<<j)

def row_bound(a,b,v,p,g):
    curves=[]
    num=(2*b+v)*v*g;den=2*b*p
    for k in (3,6,12,24,48):
        n=2*k+1;alpha=det_factor(n);base=pow_lower(g,k)
        denom2=alpha*a*base
        curves.append(ceildiv((den+num)*(n*denom2+2*k*v*b*(100<<16)),den*denom2))
    return curves

def bill(jwidth,kwidth,area,b,p):
    D=(max(1,area).bit_length()+4)//3
    best=None
    for j in range(max(20,D-3),D+3):
        if (1<<(j+1))>=p:continue
        t=prime_count(j);lost=(b.bit_length()-1)//j
        if t<=lost:continue
        c=4*(t*area+t*t*((1<<(2*j+2))+jwidth+kwidth))
        val=ceildiv(c,(t-lost)**2)
        if best is None or val<best:best=val
    z=all_moduli(jwidth,kwidth,area,b,p)
    if z is not None:best=z if best is None else min(best,z)
    extra=optimized(jwidth,kwidth,area,b,p)
    if best is None:return extra
    return min(best,extra) if extra is not None else best

@lru_cache(maxsize=20000)
def cap(a,b,v,pmin,pmax):
    if b<=a or v<=0:return v+1
    T=max(1,isqrt(pmax//b));assert T<pmin
    G=max(b*T,ceildiv(pmax,T));base=row_bound(a,b,v,pmin,b)
    best=v+1
    c=F((2*b+v)*v,2*b*pmin)
    for s in (2,4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,80):
        g=max(b,G>>s)
        lower=min(max(x,y).__ceil__() for x,y in zip(base,row_bound(a,b,v,pmin,g))) if g>b else 0
        if lower>=best:continue
        jj=c*G; kmax=F(v*b,g)+jj/2
        aw=c*v*b+jj*jj/2+3*(jj+kmax)+9
        def sieved(jmult,kmult,kj):
            U=jj*jmult;V=F(v*b,g)*kmult+jj*kj
            A=jmult*kmult*c*v*b+jmult*kj*jj*jj
            caps=[]
            for W,count in ((1,1),(3,3)):
                if W==3 and a<=3:continue
                x=bill((U/W+3).__ceil__(),(V/W+3).__ceil__(),(A/(W*W)+3*(U+V)/W+9).__ceil__(),b,pmin)
                if x is not None:caps.append(count*x)
            return min(caps) if caps else None
        case1=sieved(F(1,2),F(1),F(1,2));case2=sieved(F(1),F(1,2),F(1,2));part=max(case1,case2) if None not in (case1,case2) else None
        if part is not None:best=min(best,max(lower,part))
    return best
