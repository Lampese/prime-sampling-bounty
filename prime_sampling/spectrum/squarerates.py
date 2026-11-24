import ctypes as ct
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from fractions import Fraction as Q
from math import isqrt

from ..rectangles.old_fibre import cap as previous_cap
from ..rectangles.allsquarecap import cap as square_cap
from .fullprice import root,dirichlet_cap,ceildiv,WHEELS,weighted_near

_here=Path(__file__).resolve().parent
_tmp=TemporaryDirectory()
_bin=Path(_tmp.name)/'choices.so'
subprocess.run([os.getenv('CXX','c++'),'-O3','-std=c++17','-shared','-fPIC',str(_here/'choice.cpp'),'-o',str(_bin)],check=True)
_lib=ct.CDLL(str(_bin))
_call=_lib.choose
_call.argtypes=[ct.c_double]*10+[ct.POINTER(ct.c_int)]
_call.restype=None
sizes=((1,64),(1,16),(1,4),(1,2),(1,1),(2,1),(4,1),(8,1))

def rates(a,b,v,s,k,span,pm,p):
    sq=isqrt(p)+1; d0=(span+1)//2+1
    c=min((pm-1)//2,max(2*d0,isqrt(max(1,2*sq*d0*s))))
    best=min(Q(1),root(Q(c,c-d0)*(Q(c,s)+2*sq)/k)) if c>d0 and 2*c<pm else Q(1)
    if best<Q(1,3000):return best,best
    md=min(dirichlet_cap(a,b,v,pm,p),previous_cap(v,b,pm),square_cap(a,b,v,pm,p))
    ids=(ct.c_int*4)()
    _call(*map(float,(a,b,v,s,k,span,pm,p,md,best)),ids)
    result=[best,best]
    for code in set(ids)-{-1}:
        assert 0<=code<128
        wi,rem=divmod(code,16);mode,ix=divmod(rem,8)
        W,phi=WHEELS[wi];num,den=sizes[ix]
        D=ceildiv(span,W)+3;M=max(1,num*D//den);H=D+M
        support=M*(phi if mode else 1)
        if support<10**6 or 2*W*H>=pm:continue
        assert a>W
        J=H if mode else phi*H
        K=H if mode else W*H
        den=max(1,isqrt(v*pm//(2*b*K)));assert den<pm
        size=b*(ceildiv(p,den)+den*K);power=a*a;factors=1
        while power<=size:factors+=1;power*=a
        near=factors*(1+ceildiv(v,den)+ceildiv(2*b*den*K,pm))
        if v<=b:near=min(near,weighted_near(a,b,v,pm,p,1 if mode else W,H))
        co=(1+Q(md*near,s))/(J*s*s);rr=Q(1,support);front=Q(phi*H,k)
        const=min(Q(30001,10000)*p, p*(2*sq+1)*rr+4*sq+5+6*rr*(2*sq+6)**2)
        for j in (0,1):
            var=const+(j+1)*(p*rr)**2+2*p*rr*(1+2*rr)
            result[j]=min(result[j],root(front*(Q(1,s)+root(rr+root(co*var)))))
    return result
