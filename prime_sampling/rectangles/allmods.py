from functools import lru_cache
from math import prod,isqrt
from fractions import Fraction as F

from ..outward import b,eup
from .lattice import ceildiv
P=[q for q in range(2,10001) if all(q%i for i in range(2,isqrt(q)+1))]
SCALE=10**12
@lru_cache(None)
def local(j,y):
    D=1<<j;head=[p for p in P if p<=y];other=[p for p in P if p>y];rho=prod(F(p-1,p)for p in head);err=2**len(head)
    s2=F(sum(ceildiv(SCALE,p*p)for p in other),SCALE)+F(1,10000)
    s3=F(sum(ceildiv(SCALE**2,p*isqrt(p*SCALE**2))for p in other),SCALE)+F(2,100)
    harm=F(sum(ceildiv(SCALE,p)for p in other),SCALE)
    ld=b(F(D),'ln',True);den=b(F(10000),'ln',False)-F('1.1')
    assert D>=65536
    den0=b(F(10000),'ln',False)
    harm+=F('1.25506')*(1/b(F(D),'ln',False)+b(ld/den0,'ln',True))
    rtd=isqrt(D)+1
    low=rho*D/2*(1-2*s2)-err*(1+rtd+harm)-1
    if low<=0:return None
    D2=(rho*D/3*eup(3*s2)+err*eup(3*harm)).__ceil__()
    M=(rho*D**2/6*eup(s3)+err*D*rtd).__ceil__()
    upper_c=(rho*D/2*eup(s3)+2*err*rtd).__ceil__()
    gcd_cost=upper_c**2
    return low.__floor__(),D2,M,gcd_cost

@lru_cache(None)
def candidates(bit,top):
    out=[]
    for e in range(max(20,bit-1),min(top,bit+5)):
        for y in (3,5,7,11,13,17,19):
            row=local(e,y)
            if row is None:continue
            low,d2,ms,gc=row
            for q in (2,3,4,6):
                coef=isqrt((q+1)**3*10**4)//q+1
                a=100*(q+1)**2*d2;b=100*q*q*(q+1)*ms*ms;c=q*q*coef*gc;d=100*q*q*low**2
                out.append((a/d,b/d,c/d,a,b,c,d))
    return out

def upper(j,k,area,B,p):
    bit=(max(1,area).bit_length()+2)//3
    opts=candidates(bit,(p-1).bit_length());width=j+k
    if not opts:return None
    x=float(area);w=float(width)
    vals=[a*x+b+c*w for a,b,c,*_ in opts]
    t=min(vals)*(1+1e-13)
    return min(ceildiv(a*area+b+c*width,d) for val,(*_,a,b,c,d) in zip(vals,opts) if val<=t)
