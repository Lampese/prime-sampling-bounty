from fractions import Fraction as F
from math import prod
import subprocess

from . import certificate as c
from .mixed import col, penalties, primes

LAM=F(1,32); assert 0<LAM<=F(1,2)
def mixed_all(exe,slabs,y,mesh,rho_den):
 unit=10**16; levels=mesh*360+1;rho=F(1,rho_den)
 ps=primes(int(c.eu(y))); assert int(c.eu(y))==int(c.trans(y,'exp',False))
 norm=4+c.lu(2/rho); arr=[]
 for l,h,R in slabs:
  inv=(l+1)*c.eu(-l)
  t=(inv+2*norm*(l+1)*c.eu(-l/2))/(1-LAM-rho)
  assert 1/rho+2<2*c.trans(l/2,'exp',False)
  arr.append((col(t,mesh,levels),col(inv/(1-LAM),mesh,levels)))
 pen=' '.join(map(str,penalties(mesh,levels,unit)))
 def rec(pl,cols):
  cols=sorted(set(cols)); pp=[(p,(p**mesh-1).bit_length())for p in pl]
  dat=[f'{levels} {len(pl)} {len(cols)} {unit}',pen]
  dat += [f'{len(pl)} {j}'for j in cols]+[f'{p} {step}'for p,step in pp]
  out=subprocess.check_output([exe],input='\n'.join(dat)+'\n',text=True)
  got={}
  for row in out.splitlines():
   n,j,low,high=map(int,row.split());assert n==len(pl)
   got[j]=F(low,unit),F(high,unit)
  return got
 sc=rec(ps,[i for _,i in arr])
 ds={}
 for r in sorted(set(r for _,_,r in slabs)):
  mix=sorted(ps+c.ps[:r])
  ds[r]=rec(mix,[arr[i][0]for i,t in enumerate(slabs)if t[2]==r])
 args=[]
 for i,(l,h,r) in enumerate(slabs):
  z=h-c.ld(h-F('1.1'));g=prod(1-F(1,p)for p in c.ps[:r])
  cost=2*z*((1-LAM)*sc[arr[i][1]][1]-(1-LAM-rho)*ds[r][arr[i][0]][0]/g);assert cost>=0
  args.append((l,h,r,cost,y))
 return args
