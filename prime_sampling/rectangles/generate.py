from fractions import Fraction as F
from functools import lru_cache

from ..spectrum.modes import install_mode


def emit(left, right, mode, grid):
    from . import boxes as z
    from .sieved_lift import install as install_collision_bounds

    z.initialize(left, right)
    install_collision_bounds(z)
    install_mode(mode, z)
    base_joint = z.joint

    @lru_cache(None)
    def complete_rectangle(u,v,cl,et):
        rg,rq,mass = base_joint(u,v,cl,et)
        if not mass or u >= F(3,8) or rg < F(1,5000):
            return rg,rq,mass
        eps,_,a,b,s,width,ylo,yhi = z.params(u,v)
        eps = max(eps,z.thetaerr((1-u)*z.zl+z.b(cl,'ln',False)))
        nb = z.abel.lower_count(cl,et,eps,(1-u)*z.zh,z)
        k = (nb*ylo).__ceil__()
        assert min(k,s) >= 2
        x,y = (cl*ylo).__floor__(),(et*yhi).__ceil__()
        span = ((et-cl)*yhi).__ceil__()+2
        g1,q1 = z.rates(a,b,width,s,k,span,z.pmin,z.pmax)
        g2,q2 = z.rates(x,y,min(y-x,span),k,s,width,z.pmin,z.pmax)
        return min(rg,g1,g2),min(rq,q1,q2),mass

    z.joint = complete_rectangle
    scale = 10**12
    for index in range(67*grid//200,grid//2):
        u,v = F(index,grid),F(index+1,grid)
        cl,et = z.inner(u,v)
        rg,rq,lower = z.joint(u,v,cl,et)
        og,oq,excess = z.outer(u,v)
        masses = (z.mass(u,v),z.collar(u,v,cl,et),excess)
        print(index,*map(z.upint,(rg,rq,og,oq)),
              *((scale*x).__ceil__() for x in masses),(scale*lower).__floor__())
