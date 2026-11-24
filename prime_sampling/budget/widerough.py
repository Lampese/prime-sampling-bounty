from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, sqrt

from . import certificate as c, direct
from .certificate import dn, eta, eu, ld, lu, root, trans, up
from .direct import stride, switch

LAM=F(1,32); ALOG=c.lu(1/LAM); assert 0<LAM<=F(1,2)
@lru_cache(None)
def selweights(K):
 phi=list(range(K+1));sq=bytearray([1])*(K+1);sig=[1]*(K+1)
 for p in range(2,K+1):
  if phi[p]==p:
   for j in range(p,K+1,p):phi[j]-=phi[j]//p;sig[j]*=1+p
   sq[p*p::p*p]=bytes((K-p*p)//(p*p)+1) if p*p<=K else b''
 unit=10**12
 gd=sum(unit//phi[j] for j in range(1,K+1) if sq[j]);hu=sum((unit*sig[j]+phi[j]-1)//phi[j] for j in range(1,K+1) if sq[j])
 return F(gd,unit),F(hu,gd)

@lru_cache(None)
def weight(k):
   if k<=1:K=1;G=F(1);SM=F(1)
   else:
    K=int(trans(k,'exp',False));G=ld(K+1);SM=F(K)
    if k>=13:
     assert ld(K)+F('1.3325')-F('3.95')/isqrt(K)>k+F('1.32')
     G=max(G,k+F('1.32'))
     SM=min(SM,(F('1.53')*K+85*root(K))/G)
    if k<13 and K<=450000:
     gx,sx=selweights(K)
     G=max(G,gx);SM=min(SM,sx)
   return K,G,SM
def eta_for(u,base):
    for threshold,error in ((60,'1.232e-11'),(50,'1.135e-9'),(40,'2.149e-8'),(35,'4.977e-8'),(30,'6.076e-7'),(25,'7.437e-6'),(20,'9.164e-5')):
        if u>=threshold:return min(F(error),base)
    return base
def geometry(u,v):
    count,harm=direct.geometry(u,v)
    eta=eta_for(u,F(1))
    dx=c.eu(v-u)
    return min(count,c.up((dx-1+(dx+1)*eta)/u)),min(harm,c.up((1+eta)*c.lu(v/u)+2*eta/u))
def bound(l,h,R,ylog,rough,dw,alphas,*,options_for_rank):
 zl=dn(l-lu(l+1));zh=up(h-ld(h-F('1.1')))
 invx=(l+1)*eu(-l)
 px=(l+1)*eu(-l/2)
 delta,opts,_=options_for_rank(R)
 fo=[tuple(map(float,row[:3])) for row in opts]
 def cells(start,end):
  while start<end:
   v=min(end,start+(stride if start<switch else dw),switch if start<switch else end)
   yield start,v
   start=v
 def band(u,v,U,high):
  c,H=geometry(u,v);da=v-u;cb=F('.6')/U if 2*U<=zl-v else max(F('.6')/U,1/(zl-v-F('1.1'))+eu(v)*invx)
  W=up(2*zh*cb*H)
  assert zl-v>23
  if 4*U>=zh-u and U>=23:
   def hh(maxlog):return (1+eta)*lu(maxlog/U)+2*eta/U if maxlog>U else F(0)
   M=zh-u; ylo=zl-u
   ll=ALOG+da
   assert ylo/3-ll>=23
   den=1/(ylo-ll)+hh(M/2)/(ylo/2-ll)
   if 3*U<=M:den+=hh(M/3)*hh((M-U)/2)/(ylo/3-ll)
   W=min(W,up(2*zh*H*(1-LAM+(1+LAM)*eta)*den))
  else:den=None
  tilt=eu(u); coarse=cb
  diag=4*zh**2*coarse*H*eu(-u)*(1-LAM)
  off=4*zh**2*coarse*c*c*tilt*3*(h+2)*px
  params=[]
  for logk in (0,2,4,6,8,10,12,14,16,18,20,22):
   k=min(F(logk),U-F('.001'))
   K,G,SM=weight(k)
   di=diag/G+8*zh**2*coarse*c*SM**2*invx
   params.append((W,F(0),di,off*SM**2))
  if high:
   ka=min(U-F('.001'),(zl-v)/3-F('.001'))
   edge=2*zh*H*(eu(da)-1)/ka+4*zh*c*eu(u+2*ka)*invx
   if den:
    tail_eta=eta_for(min(U,(zl-u)/3-ALOG-da),eta)
    edge=min(edge,2*zh*H*(eu(da)-1+4*tail_eta*eu(da))*den)
    cb=min(cb,(1-LAM*trans(-da,'exp',False)+2*eta)*den)
   for m in (F(1,8),F(1,2),F(2),F(12),F(80)):
    length=1-LAM*trans(-da,'exp',False)+m+F('1e-8')
    assert eu(u)*invx < F('1e-9')
    af=root(length/m)+F('1e-6')
    di=4*zh**2*cb*c*eu(-u)*length
    of=4*zh**2*cb*c*c*tilt*3*af*px
    params.append((W+edge,edge,di,of))
  best=(float('inf'),0,0)
  for j,row in enumerate(params):
   ww,edge,dd,oo=map(float,row)
   for k,(gap,var,norm) in enumerate(fo):
    a=gap*ww+edge+sqrt(var*dd+norm*norm*oo)
    if a<best[0]:best=(a,j,k)
  ww,edge,dd,oo=params[best[1]];gap,var,norm,_,_=opts[best[2]]
  return up(gap*ww+edge+root(var*dd+norm**2*oo))
 return band
