from fractions import Fraction as Q
from functools import lru_cache
import sys

from ..outward import elo,eup,b,upint,nthroot
from ..spectrum.squarerates import rates
from . import abel
LAM=Q(1,32)
ALOG=b(1/LAM, "ln", True)
def initialize(left, right):
    global lo,hi,zl,zh,pmin,pmax
    lo,hi=left,right
    zl=b(Q((elo(lo)/(lo-1)).__floor__()), 'ln', False)
    zh=hi-b(hi-Q('1.1'),'ln',False)
    pmin=elo(lo).__floor__()
    pmax=eup(hi).__ceil__()
def thetaerr(z):
    assert z>=40
    return Q(2149,10**11)

@lru_cache(None)
def Xlo(t):return elo(t*(zl if t>=0 else zh))
@lru_cache(None)
def Xhi(t):return eup(t*(zh if t>=0 else zl))
@lru_cache(None)
def params(u,v):
    eps=thetaerr(u*zl)
    na=((1-eps)*Xlo(v-u)-(1+eps))/(v*zh)
    aa=Xlo(u).__floor__();bb=Xhi(v).__ceil__()
    ssize=(na*aa).__ceil__();width=(bb*(1-Xlo(u-v))).__ceil__()
    return eps,na,aa,bb,ssize,width,Xlo(1-u),Xhi(1-u)
@lru_cache(None)
def joint(u,v,cl,et):
    if cl>=et:return Q(1),Q(1),Q(0)
    eps,na,aa,bb,ssize,width,ylo,yhi=params(u,v)
    eps=max(eps,thetaerr((1-u)*zl+b(cl,'ln',False)))
    nb=abel.lower_count(cl,et,eps,(1-u)*zh,sys.modules[__name__])
    kk=(nb*ylo).__ceil__()
    if min(kk,ssize)<2 or na<=0:return Q(1),Q(1),Q(0)
    xx=(cl*ylo).__floor__();yy=(et*yhi).__ceil__();span=((et-cl)*yhi).__ceil__()+2
    sp=eup((u-Q(1,2))*lo+(1-u)*b(hi,'ln',True))
    old=min(Q(1),nthroot(Q(100,99)/nb*(101*Xhi(-u)/na+2*sp)))
    if u<Q('.375'):
        return old,old,2*zl*na*nb
    r1=rates(aa,bb,width,ssize,kk,span,pmin,pmax)
    r2=rates(xx,yy,min(yy-xx,span),kk,ssize,width,pmin,pmax)
    return min(old,r1[0],r2[0]),min(old,r1[1],r2[1]),2*zl*na*nb
def inner(u,v):
    return (max(LAM,Xhi(u+v-1)),Xlo(u-v))
def mass(u,v):
    eps=thetaerr(u*zl)
    hh=(1+eps)*b(v/u,'ln',True)+2*eps/(u*zl)
    denom=max(u, 1-v-ALOG/zl)
    tri=min(2*(1-LAM+(1+LAM)*eps),2*(1-Xlo(2*u-1)+2*eps)+Q(2,10**15))
    return min(hh*tri/denom,abel.upper_mass(u,v,sys.modules[__name__]))

def collar(u,v,cl,et,steps=8):
    eps=thetaerr(u*zl)
    if cl>=et:return mass(u,v)
    ans=Q(0);den=max(u, 1-v-ALOG/zl)
    for j in range(steps):
        left=u+(v-u)*j/steps;right=u+(v-u)*(j+1)/steps
        top=Xhi(u-left); bot=max(LAM*Xlo(u-right),Xlo(u+left-1))
        delta=max(0,min(cl,top)-bot)+max(0,top-max(et,bot))+4*eps*top
        np=Xhi(left-u)*((1+eps)*Xhi(right-left)-(1-eps))/(left*zl)
        ans+=2*np*delta/den
    return min(ans,mass(u,v))

def outer(u,v,steps=8):
    cl=max(LAM*Xlo(u-v),Xlo(2*u-1));et=Q(1)
    rg,rq,mm=joint(u,v,cl,et)
    den=1-u+b(cl,'ln',False)/zl
    eps=thetaerr(min(u*zl,den*zl))
    extra=Q(0)
    for j in range(steps):
        l=u+(v-u)*j/steps;r=u+(v-u)*(j+1)/steps
        floor_q=max(LAM*Xhi(u-l),Xhi(u+r-1))
        delta=max(0,min(1,floor_q)-cl) + max(0,1-max(cl, Xlo(u-r)))+4*eps
        np=Xhi(l-u)*((1+eps)*Xhi(r-l)-(1-eps))/(l*zl)
        extra+=2*np*delta/den
    return rg,rq,extra
