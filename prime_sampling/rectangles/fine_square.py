from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

from ..outward import b
from .lattice import ceildiv

@lru_cache(None)
def weights(D):
    sums=[0,0,0]
    low=D//8
    for i in range(1,8):
        h=D*(i+1)//8;w=8-i
        ll=b(F(low),'ln',False);hh=b(F(h),'ln',True)
        if ll>=23:eta=F('.01')/ll**2
        elif ll>=19:eta=F('.05')/ll**2
        elif ll>=16:eta=F('.2')/ll**2
        else:eta=None
        if eta is not None:
            lower=(((1-eta)*h-(1+eta)*low)/hh).__ceil__()
            upper=(((1+eta)*h-(1-eta)*low)/ll).__floor__()
        else:
            lower=(h/(hh-1)-low/(ll-F('1.1'))).__ceil__()
            upper=(h/(b(F(h),'ln',False)-F('1.1'))-low/(b(F(low),'ln',True)-1)).__floor__()
        assert upper>=lower>=1
        sums[0]+=w*lower
        sums[1]+=w*w*upper
        sums[2]+=w*h*upper
        low=h
    return sums

def optimized(jwidth,kwidth,area,B,p):
    D=(max(1,area).bit_length()+3)//3
    best=None
    for j in range(max(23,D-1),D+6):
        T=1<<j
        if T>=p:continue
        t,w2,mp=weights(T)
        bad=7*((B.bit_length()-1)//(j-3))
        if t<=bad:continue
        upper_t=weights_total_up(T)
        val=ceildiv(4*(w2*area+mp*mp+upper_t**2*(jwidth+kwidth)),(t-bad)**2)
        for q in (2,4,6,8):
            if max(jwidth,kwidth)*(q+1)>=q*(T//8)**2:continue
            ans=ceildiv((q+1)**2*area*w2+mp**2*(q+1)*q*q, q*q*(t-bad)**2)
            val=min(val,ans)
        if best is None or val<best:best=val
    return best

@lru_cache(None)
def weights_total_up(D):
    upper=(F(D)/(b(F(D),'ln',False)-F('1.1'))).__floor__()
    return isqrt(upper*weights(D)[1])+1
