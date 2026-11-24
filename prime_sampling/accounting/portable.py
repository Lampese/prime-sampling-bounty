from fractions import Fraction as F
from functools import cache
from math import isqrt

from ..outward import b, elo, eup

UNIT=10**12
LAMBDA=F(1,32)
K=F(11,10)

@cache
def scale(left,right):
    assert 40<=left<right<=175
    x0=(elo(left)/(left-1)).__floor__()
    low=b(F(x0),'ln',False)
    high=right-b(right-K,'ln',False)
    assert min(elo(low)/32,elo(F(67,200)*low),elo(low/2))>60184
    return x0,low,high

def prime_floor(left,right):
    _,low,high=scale(left,right)
    k=b(F(32),'ln',True)+K
    bound=2*high/(high-1)-2*LAMBDA*low/(low-k)
    return (UNIT*bound).__floor__()

def rough_ceiling(left,right):
    x0,low,high=scale(left,right)
    primes=[n for n in range(2,200) if all(n%d for d in range(2,isqrt(n)+1))]
    g=F(1); opts=[1-LAMBDA+F(1,x0)]
    for j,p in enumerate(primes,1):
        assert p<elo(F(67,200)*low)
        g*=1-F(1,p)
        opts.append((1-LAMBDA)*g+F(1<<j,x0))
    return (2*high*min(opts)*10**9).__ceil__()

def enclosing_boxes(left,right,path,grid=600):
    _,low,high=scale(left,right)
    assert grid>=600 and grid%200==0
    with path.open('w') as dst:
        for i in range(67*grid//200,grid//2):
            u,v=F(i,grid),F(i+1,grid)
            side=eup((v-u)*high)/(v*low-K)-1/(u*high-1)
            mass=(2*high*max(F(0),side)/((1-u)*low-K)*UNIT).__ceil__()
            print(i,*([10**10]*4),mass,mass,mass,0,file=dst)
