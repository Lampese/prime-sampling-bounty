from fractions import Fraction as F
from math import prod,isqrt
from pathlib import Path
from itertools import accumulate
import ctypes as ct,subprocess,tempfile,os
SOURCE=Path(__file__).resolve().with_suffix('.cpp')
def nroot(x,k):
 unit=10**20
 n=(x*unit**k).__ceil__()
 if n<=1:return F(n,unit)
 y=1<<((n.bit_length()+k-1)//k)
 while True:
  nxt=((k-1)*y+n//y**(k-1))//k
  if nxt>=y:return F(y+1,unit)
  y=nxt
def install(z):
 from . import squarerates as sr
 tmp=tempfile.TemporaryDirectory()
 lib=Path(tmp.name)/'split.so'
 subprocess.run([os.getenv('CXX','c++'),'-O3','-std=c++17','-shared','-fPIC',SOURCE,'-o',lib],check=True)
 linked=ct.CDLL(str(lib)); select=linked.choose;select.argtypes=[ct.c_double]*10+[ct.POINTER(ct.c_int)];select.restype=None
 primes=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47)
 rows=[(prod(primes[:d]),prod(p-1 for p in primes[:d]),prod(primes[:i]),prod(p-1 for p in primes[:i]),2**(d-i))for d in range(1,len(primes)+1)for i in range(min(d,5)+1)]
 ratios=[F(x) for x in ('1/128','1/32','1/8','1/3','.7','1.5','3','6','12','24')]
 def rates(a,b,v,s,k,span,pm,p):
  orig=zsaved(a,b,v,s,k,span,pm,p)
  if max(orig)<F(1,3000):return orig
  md=min(sr.dirichlet_cap(a,b,v,pm,p),sr.previous_cap(v,b,pm),sr.square_cap(a,b,v,pm,p))
  res=list(orig);ids=(ct.c_int*2)();select(*map(float,(a,b,v,s,k,span,pm,p,md,max(orig))),ids)
  root=sr.root;sqrt=isqrt(p)+1
  for cid in set(ids)-{-1}:
   i,j=divmod(cid,10);W,phi,W1,p1,divs=rows[i];p2=phi//p1
   assert W<a and W%W1==0
   D=sr.ceildiv(span,W)+5;M=max(1,(D*ratios[j]).__floor__());H=D+2*M;S=M*p2;J=H*p1
   if S<10**6 or 2*W*H>=pm:continue
   near=sr.weighted_near(a,b,v,pm,p,W1,H)
   co=(1+F(md*near,s))/(J*s*s);rr=F(2*M*M+1,3*S*M*M);front=F(phi*H,k)
   l1sq=min(F(1,p2),(F(divs,p2)*(1+root(F(2*(W//W1),M))))**2)
   c=14*F(p*p*sqrt,M*M)*rr*l1sq
   for q in (0,1):
    center=rr if q else rr/2
    var=(c+2*(p*rr)**2+2*p*rr*(1+2*rr)) if q else (c+(p*rr)**2/4+8*p*rr*(1+rr))
    res[q]=min(res[q],root(front*(F(1,s)+root(center+root(co*var)))))
    for r,mult,factorial in ((3,15,6),(4,105,24),(5,945,120)):
     total=mult*(2**r if q else F(1,2**r));lead=mult if q else F(factorial,2**r)
     mixed=(4**r-2)*p*sqrt*rr*(F(p,M*M)*l1sq)**(r-1) +lead*p*p*rr**r+(total-lead)*p*rr**(r-1)
     res[q]=min(res[q],root(front*(F(1,s)+nroot(co*mixed,2*r))))
  return res
 zsaved=z.rates;z.rates=rates
