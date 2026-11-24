from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess

from . import band, menu, tight, wideseed
from . import widerough as s

HERE = Path(__file__).resolve().parent

def build(left,right,rank,method,size,fine):
    if str(rank) not in tight.patterns:
        tight.patterns[str(rank)]=[(h,n,t) for h in range(min(rank,16)+1)
                                  for n in range(h,rank+1) for t in (0,2,4) if t<=n-h]
    density,options,_=menu.group(rank)
    prices=[row[:3] for row in options]
    choices=prices,[tuple(map(float,row)) for row in prices]
    z=s.up(right-s.ld(right-F('1.1')))
    transfer=s.bound(left,right,rank,F(13),0,fine,[],options_for_rank=menu.group)
    seeds=((F(12),40,32768),(F('12.5'),48,65536),(F(13),48,65536),
           (F('13.5'),36,65536),(F(14),48,65536)) if method=='multi' else ((F(13),24,8192),)
    with TemporaryDirectory() as tmp:
        program=Path(tmp)/'rec'
        subprocess.run([os.getenv('CXX','c++'),'-O3','-std=c++17',HERE/'rec.cpp','-o',program],check=True)
        starts={}
        for start,mesh,rho in seeds:
            value=wideseed.mixed_all(program,[(left,right,rank)],start,mesh,rho)[0][3]
            starts[start]=value+2*z*(1/density-1)/(s.trans(start,'exp',False)-1)
    if method=='multi':
        edges=[F(i,2) for i in range(24,53)]
        best={}
        for i,v in enumerate(edges):
            candidates=[starts[v]] if v in starts else []
            candidates += [best[u]+transfer(u,v,u,False) for u in edges[max(0,i-7):i]]
            best[v]=min(candidates)
        total=best[F(26)]
    else:
        u=F(13);total=starts[u]
        while u<26:
            v=min(F(26),u+F('1.5'))
            total+=transfer(u,v,u,False)
            u=v
    u=F(26);out=[]
    for cut in (67,68,69,70):
        end=z*cut/200
        while u<end:
            v=min(u+size,end)
            count=max(1,((v-u)/fine).__ceil__())
            total+=count*band.price(band.params(left,right,u,v,count),choices)
            u=v
        out.append((total*10**9).__ceil__())
    return out
