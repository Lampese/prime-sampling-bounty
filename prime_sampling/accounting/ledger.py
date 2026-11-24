from fractions import Fraction as F
from functools import cache
from math import isqrt, lcm, prod
from pathlib import Path
from tempfile import TemporaryDirectory

from ..outward import b, eup
from . import portable
from .catalog import BUDGETS, TABLES

I, J = 10**10, 10**12
PRIMES = tuple(p for p in range(2,256) if all(p%d for d in range(2,isqrt(p)+1)))
CUTS = (67,68,69,70)
TABLE_HEADER = 'strip\tg_in\tq_in\tg_out\tq_out\tmass_max\tcollar_max\texcess_max\tinside_min'

def ceil(n,d):
    return -(-n//d)

def fmt(n):
    sign = '-' if n < 0 else ''
    return f'{sign}{abs(n)//J}.{abs(n)%J:012d}'

@cache
def allowance(upper):
    assert 150 <= upper <= 175
    w,eta = F(1,32),F(2454,10**15)
    a = b(1/w,'ln',True)
    z = upper-b(upper-F('1.1'),'ln',False)
    m = 1-w-(1+w)*eta
    j = (1-eta)*(1-w)-(1+eta)*w*a
    k = 2*(1-eta)*(1-w-w*a)-(1+eta)*w*a*a
    assert min(j,k)>0
    return (2*(m+j/z+k/z**2)*J).__floor__()

@cache
def rank(upper,branch):
    seq = tuple(p for p in PRIMES if branch!='without3' or p!=3)
    bound = eup(upper)
    r = next(i for i in range(1,len(seq)) if prod(seq[:i+1])>bound)
    return seq[:r]

@cache
def bank(seq,kind):
    density = prod(1-F(1,p) for p in seq)
    scale = (J/density).__ceil__()
    entries = set()
    for h in range(0 if kind=='single' else 2 if kind=='six' else 1, min(12,len(seq))+1):
        prefix = prod(1-F(1,p) for p in seq[:h])
        E = [F(1)]+[F(0)]*4
        C = E.copy()
        inv = F(0)
        for t in range(h,len(seq)+1):
            if t>h:
                p = seq[t-1]
                inv += F(1,p)
                if inv>1:
                    break
                for i in range(4,0,-1):
                    E[i] += E[i-1]/p
                    C[i] += (1-F(2,p))*C[i-1]
            if kind=='single' and t!=len(seq) and t-h not in (0,1,2,3,4,6,8,12,16,20):
                continue
            for degree in (0,2,4):
                mean = prefix*sum((-1)**i*E[i] for i in range(degree+1))
                norm = prefix*2**h*sum(C[:degree+1])-mean
                m = ((1-mean)*I).__floor__()
                if kind=='single':
                    entries.add((m,(norm*I).__ceil__(),0))
                elif mean<=1:
                    inner = (3 if kind=='six' else 1)*mean
                    entries.add((m,((norm-inner)*I).__ceil__(),(inner*I).__ceil__()))
    frontier = []
    for m,n,q in sorted(entries,key=lambda row:(-row[0],row[1],row[2])):
        if not any(o>=m and v<=n and u<=q for o,v,u in frontier):
            frontier.append((m,n,q))
    return tuple(frontier),scale

@cache
def check_rows(file,grid):
    assert grid>=600 and grid%200==0
    lines = file.read_text().splitlines()
    if lines and lines[0] == TABLE_HEADER:
        lines = lines[1:]
    result=[]
    for j,line in enumerate(lines,67*grid//200):
        row=tuple(map(int,line.split()))
        assert len(row)==9 and row[0]==j
        assert all(0<=v<=I for v in row[1:5]) and min(row[5:])>=0
        result.append(row[1:])
    assert len(result)==33*grid//200
    return result

@cache
def strip_costs(file,grid,kind,seq):
    options,scale=bank(seq,kind)
    def excluded(g,q):
        if kind=='single':
            return max(0,max(m-ceil(n*g+(I-m)*max(0,q-g),I) for m,n,_ in options))
        return max(0,max(m-ceil(n*g+v*q,I) for m,n,v in options))
    result=[]
    for g,q,go,qo,m,h,e,m0 in check_rows(file,grid):
        fi,fo=excluded(g,q),excluded(go,qo)
        cap=min(m-(fi*m0)//I,ceil((I-fi)*m+fi*h,I),ceil((I-fo)*(m+e),I))
        assert cap>=0
        result.append(ceil(scale*cap,J)-m)
    return result

def suffix_costs(selected,seq,tiles):
    rows=[(g,strip_costs(path,g,k,seq)) for _,k,_,_,g,path in selected]
    if not tiles:
        macros=[[sum(cost[j*g//200:(j+1)*g//200]) for j in range(33)] for g,cost in rows]
        best=list(map(min,zip(*macros)))
        return {cut:sum(best[cut-67:]) for cut in CUTS}
    unit=lcm(*(g for g,_ in rows))
    edges={}
    for g,cost in rows:
        stride=unit//g
        for j,v in enumerate(cost,g*67//200):
            a=j*stride
            ends=edges.setdefault(a,{})
            ends[a+stride]=min(v,ends.get(a+stride,v))
    tail={unit//2:0}
    for a,ends in sorted(edges.items(),reverse=True):
        tail[a]=min(v+tail[end] for end,v in ends.items())
    return {cut:tail[cut*unit//200] for cut in CUTS}

def evaluate(a,z,stage,branch,extras=(),added_budgets=(),main_lower=None):
    A = allowance(z) if main_lower is None else main_lower
    seq=rank(z,branch)
    groups={'A':{'A'},'B':{'B','shared'},'C':{'shared'},'D':{'D'},
            'E':{'shared'},'probe':{'shared','extension','B','D','A','extra'}}[stage]
    candidates=[v for s,l,u,r,v in [*BUDGETS,*added_budgets] if s in groups and l<=a<z<=u and r>=len(seq)]
    sources={'A':{'A','prior'},'B':{'B'},'C':{'C'},'D':{'D','prior'},'E':{'E'},
             'probe':{'A','prior','B','C','D','E','extra'}}[stage]
    selected=[row for row in [*TABLES,*extras] if row[0] in sources and row[2]<=a<z<=row[3]
              and (row[1]!='six' or branch=='with3')
              and not (stage=='C' and branch=='with3' and row[1]!='six')]
    if not candidates or not selected:
        raise ValueError(f'No common input for [{a},{z}], {branch}. Split at an input boundary.')
    charges=suffix_costs(selected,seq,stage=='D')
    cuts=(67,68,70) if stage in ('A','D') else CUTS
    results=[]
    for cut in cuts:
        bounds=[row[cut-67] for row in candidates if row[cut-67] is not None]
        if bounds:
            R=min(bounds)*1000
            results.append((A-R-charges[cut],cut,R,charges[cut]))
    value,cut,R,T=max(results)
    return len(seq),value,cut,R,T

def report(a,z,branch,result):
    r,m,j,R,T=result
    print(f'[{a}, {z}] {branch}  r={r}  d={j}/200  A={fmt(m+R+T)}  R={fmt(R)}  T={fmt(T)}  margin={fmt(m)}')

def probe(a,z,tables=(),budget_files=(),new_lower=None):
    extras=list(tables)
    additions=[]
    for file in budget_files:
        for line in file.read_text().splitlines():
            l,u,r,*values=line.split()
            l,u,r=F(l),F(u),int(r)
            values=list(map(int,values))
            if not (40<=l<u<=175 and 1<=r<=43 and len(values)==4 and min(values)>=0):
                raise ValueError(f'Invalid budget row in {file}.')
            additions.append(('extra',l,u,r,values))
    main_lower = allowance(z) if a>=150 else portable.prime_floor(a,z)
    if new_lower is not None:
        main_lower=max(main_lower,(J*new_lower).__floor__())
    with TemporaryDirectory() as tmp:
        if a<150:
            raw=Path(tmp)/'outer.dat'
            portable.enclosing_boxes(a,z,raw)
            extras.append(('extra','single',a,z,600,raw))
            R=portable.rough_ceiling(a,z)
            additions.append(('extra',a,z,43,[R]*4))
        return {branch:evaluate(a,z,'probe',branch,extras,additions,main_lower)
                for branch in ('with3','without3')}
