from fractions import Fraction as F
from math import prod,isqrt
from pathlib import Path
from itertools import accumulate
import ctypes as ct,subprocess,tempfile,os
SOURCE=Path(__file__).resolve().with_suffix('.cpp')
def install(z):
 from . import squarerates as sr
 tmp=tempfile.TemporaryDirectory()
 lib=Path(tmp.name)/'split.so'
 subprocess.run([os.getenv('CXX','c++'),'-O3','-std=c++17','-shared','-fPIC',SOURCE,'-o',lib],check=True)
 linked=ct.CDLL(str(lib)); select=linked.choose;select.argtypes=[ct.c_double]*10+[ct.POINTER(ct.c_int)];select.restype=None
 primes=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47)
 rows=[(prod(primes[:d]),prod(p-1 for p in primes[:d]),prod(primes[:i]),prod(p-1 for p in primes[:i]))for d in range(1,len(primes)+1)for i in range(min(d,5)+1)]
 ratios=[F(x) for x in ('1/128','1/32','1/8','1/3','.7','1.5','3','6','12','24')]
 def rates(a,b,v,s,k,span,pm,p):
  orig=zsaved(a,b,v,s,k,span,pm,p)
  if max(orig)<F(1,3000):return orig
  md=min(sr.dirichlet_cap(a,b,v,pm,p),sr.previous_cap(v,b,pm),sr.square_cap(a,b,v,pm,p))
  res=list(orig);ids=(ct.c_int*2)();select(*map(float,(a,b,v,s,k,span,pm,p,md,max(orig))),ids)
  root=sr.root;sqrt=isqrt(p)+1
  for cid in set(ids)-{-1}:
   i,j=divmod(cid,10);W,phi,W1,p1=rows[i];p2=phi//p1
   assert W<a and W%W1==0
   D=sr.ceildiv(span,W)+5;M=max(1,(D*ratios[j]).__floor__());H=D+2*M;S=M*p2;J=H*p1
   if S<10**6 or 2*W*H>=pm:continue
   near=sr.weighted_near(a,b,v,pm,p,W1,H)
   co=(1+F(md*near,s))/(J*s*s);rr=F(2*M*M+1,3*S*M*M);front=F(phi*H,k)
   c=16*F(p*p*sqrt,M**3)
   for q in (0,1):
    var=c+2*(p*rr)**2+2*p*rr*(1+2*rr)
    res[q]=min(res[q],root(front*(F(1,s)+root(rr+root(co*var)))))
  return res
 zsaved=z.rates;z.rates=rates
