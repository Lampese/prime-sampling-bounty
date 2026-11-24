from fractions import Fraction as F

from . import widerough as s
def params(lo,hi,a,b,n):
    w=(b-a)/n;T=b-w
    zl=s.dn(lo-s.lu(lo+1));zh=s.up(hi-s.ld(hi-F('1.1')))
    inv=(lo+1)*s.eu(-lo);px=(lo+1)*s.eu(-lo/2);eps=s.eta_for(a,F('.01')/a**2);eta=s.eta;lam=s.LAM
    U=s.eu(w);c=s.up((U-1+(U+1)*eps)/a);H=min(c,s.up((1+eps)*s.lu(1+w/a)+2*eps/a))
    cb=F('.6')/a
    if 2*T>zl-b:cb=max(cb,1/(zl-b-F('1.1'))+s.eu(b)*inv)
    W=s.up(2*zh*cb*H);den=None
    if 4*a>=zh-a:
        def hh(M):return (1+eta)*s.lu(M/a)+2*eta/a if M>a else F(0)
        M=zh-a;y=zl-T;ll=s.ALOG+w;assert y/3-ll>=23
        den=1/(y-ll)+hh(M/2)/(y/2-ll)
        if 3*a<=M:den+=hh(M/3)*hh((M-a)/2)/(y/3-ll)
        W=min(W,s.up(2*zh*H*(1-lam+(1+lam)*eta)*den))
    diag=4*zh**2*cb*H*s.eu(-a)*(1-lam);off=4*zh**2*cb*c*c*s.eu(T)*3*(hi+2)*px
    rows=[]
    for kk in range(0,23,2):
        k=min(F(kk),a-F('.001'));K,g,sm=s.weight(k)
        di=diag/g+8*zh**2*cb*c*sm*sm*inv;rows.append((W,F(0),di,off*sm*sm))
    ka=min(a-F('.001'),(zl-b)/3-F('.001'))
    edge=2*zh*H*(U-1)/ka+4*zh*c*s.eu(T+2*ka)*inv
    if den is not None:
        ep=s.eta_for(min(a,(zl-T)/3-ll),eta)
        edge=min(edge,2*zh*H*(U-1+4*ep*U)*den)
        cb=min(cb,(1-lam*s.trans(-w,'exp',False)+2*eta)*den)
    for rad in (F(1,8),F(1,2),F(2),F(12),F(80)):
        length=1-lam*s.trans(-w,'exp',False)+rad+F('1e-8')
        assert s.eu(T)*inv<F('1e-9')
        af=s.root(length/rad)+F('1e-6')
        di=4*zh**2*cb*c*s.eu(-a)*length
        oo=4*zh**2*cb*c*c*s.eu(T)*3*af*px
        rows.append((W+edge,edge,di,oo))
    return rows

def price(rows,opts):
    best=(float('inf'),0,0)
    for j,(W,e,di,of) in enumerate(rows):
        w,e0,d,x=map(float,(W,e,di,of))
        for k,(gap,var,norm) in enumerate(opts[1]):
            v=gap*w+e0+(var*d+norm*norm*x)**.5
            if v<best[0]:best=(v,j,k)
    _,i,k=best;W,e,di,of=rows[i];gap,var,norm=opts[0][k]
    return s.up(gap*W+e+s.root(var*di+norm*norm*of))
