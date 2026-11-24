from fractions import Fraction as F

from . import certificate as c
switch=26
dw=F('.005')
stride=F('1.5')
def geometry(u,v):
 d=v-u;ct=c.eu(d)/(v-F('1.1'))-1/(u-1)
 hh=min(ct,1/(v-F('1.1'))-1/(u-1)+c.lu((v-F('1.1'))/(u-F('1.1'))))
 if u>=F(12):
  eta=0
  if u>=F('23'):eta=F('.01')/u**2
  elif u>=F('19'):eta=F('.05')/u**2
  elif u>=F('16'):eta=F('.2')/u**2
  if eta:
   ct=min(ct,(c.eu(d)-1+eta*(c.eu(d)+1))/u)
   hh=min(hh,(1+eta)*c.lu(v/u)+2*eta/u)
 return c.up(ct),c.up(hh)
