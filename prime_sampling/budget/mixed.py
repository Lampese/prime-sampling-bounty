from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

@lru_cache(None)
def primes(n):
 v=bytearray([1])*(n+1)
 for p in range(2,isqrt(n)+1):
  if v[p]:v[p*p:n+1:p]=bytes((n-p*p)//p+1)
 return [p for p in range(2,n+1) if v[p]]
def col(x,mesh,levels):
 a=1/x**mesh
 j=a.numerator.bit_length()-a.denominator.bit_length()
 if j>0 and F(2)**j>a:j-=1
 return min(levels-1,max(0,j))
@lru_cache(None)
def penalties(mesh,levels,unit):
 def penalty(j):
  lo,hi=0,unit
  while hi-lo>1:
   m=(hi+lo)//2
   if (m**mesh<<j)>=unit**mesh:hi=m
   else:lo=m
  return hi
 return [penalty(j)for j in range(levels)]
