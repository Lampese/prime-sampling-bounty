from functools import lru_cache
from fractions import Fraction as F
from math import isqrt

def ceildiv(a,b): return -(-a//b)

@lru_cache(None)
def grid_power(j,k):
    n=2*k+1;target=1<<(j*k+16*n)
    lo=1<<((target.bit_length()-1)//n);hi=2*lo
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**n<=target:lo=mid
        else:hi=mid
    return lo

@lru_cache(None)
def det_factor(n):
    return next(t for t in range(400,99,-1) if t**(n*(n-1))*n**n*2**(2*(n-1))<=400**(n*(n-1)))

@lru_cache(maxsize=100000)
def dirichlet_cap(a,b,v,pmin,pmax):
    if b*b>25*pmax or b<=a or v<=0:return v+1
    T=max(1,isqrt(pmax//b));assert T<pmin
    cap=max(T*b,ceildiv(pmax,T))+ceildiv(T*v,2)
    assert cap>=b
    best=v+1
    for k in (3,6,12,24,48):
        n=2*k+1
        factor=det_factor(n)
        upper=0
        for G in (b,cap):
            j=G.bit_length()-1
            powlo=(grid_power(j,k)*((2<<j)-G)+grid_power(j+1,k)*(G-(1<<j)))//(1<<j)
            assert powlo>0
            per=F(n)+F(2*k*v*b*(100<<16),factor*a*powlo)
            upper=max(upper,((1+F(G*v,pmin))*per).__ceil__())
        best=min(best,upper)
    return best
