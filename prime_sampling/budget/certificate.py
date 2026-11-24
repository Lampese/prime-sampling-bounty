from fractions import Fraction as F
from functools import lru_cache
from decimal import Decimal,Context,ROUND_CEILING,ROUND_FLOOR
from math import comb,prod,isqrt
Q=10**30; eta=F(1,50000); zmin=13; dw=F(3,200); switch=26
def up(a):return F((a*Q).__ceil__(),Q)
def dn(a):return F((a*Q).__floor__(),Q)
def root(a):return F(isqrt((a*Q*Q).__ceil__())+1,Q) if a else F(0)
@lru_cache(None)
def trans(a,op,upper):
 a=F(a)
 if a==0 and op=='exp':return F(1)
 if a==1 and op=='ln':return F(0)
 c=Context(prec=51,rounding=ROUND_CEILING if upper else ROUND_FLOOR)
 v=c.divide(Decimal(a.numerator),Decimal(a.denominator))
 vv=getattr(c,op)(v)
 return F(c.next_plus(vv) if upper else c.next_minus(vv))
def eu(a):return trans(a,'exp',True)
def ld(a):return trans(a,'ln',False)
def lu(a):return trans(a,'ln',True)
ps=[p for p in range(2,251) if all(p%k for k in range(2,isqrt(p)+1))]
@lru_cache(None)
def group(R):
 delta=prod(1-F(1,p) for p in ps[:R]); upper=[]; lower=[]
 for h in range(0,R+1):
  head=prod(1-F(1,p) for p in ps[:h]);tail=ps[h:R];dist=[F(1)]
  for p in tail:
   q=F(1,p)
   dist=[a*(1-q)+b*q for a,b in zip(dist+[0],[0]+dist)]
  for t in (0,2,4,6,8,len(tail) if len(tail)%2==0 else len(tail)+1):
   vals=[1 if n==0 else (comb(n-1,t) if n>t else 0) for n in range(len(dist))]
   mean=head*sum(prob*k for prob,k in zip(dist,vals))/delta
   var=head*sum(prob*k*k for prob,k in zip(dist,vals))/delta**2-mean*mean
   norm=2**h*head/delta*sum(comb(len(tail),j) for j in range(min(t,len(tail))+1))
   upper.append((up(mean-1),up(var),up(norm),h,t))
   if t not in (0,2,4,6,8):continue
   j=t+1
   mu1=mean
   mu2=head*sum(prob*comb(n,j) for n,prob in enumerate(dist) if n>=j)/delta
   norm2=2**h*head/delta*(comb(len(tail),j) if j<=len(tail) else 0)
   lower.append((mu1,mu2,norm,norm2,h,j))
 vals=[]
 for b in sorted(set(upper)):
  if not any(a[0]<=b[0] and a[1]<=b[1] and a[2]<=b[2] for a in vals):vals.append(b)
 return delta,vals,lower
