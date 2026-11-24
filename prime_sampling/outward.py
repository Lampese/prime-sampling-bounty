from fractions import Fraction as Q
from math import isqrt,prod
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR
from functools import lru_cache
@lru_cache(None)
def b(x,f,side):
    c=Context(prec=40,rounding=ROUND_CEILING if side else ROUND_FLOOR)
    t=c.divide(Decimal(x.numerator),Decimal(x.denominator))
    y=getattr(c,f)(t)
    return Q(c.next_plus(y) if side else c.next_minus(y))

def elo(x): return b(x, 'exp', False)
def eup(x): return b(x, 'exp', True)
def upint(x): return (x*10**10).__ceil__()
def nthroot(x): return Q(isqrt((x*10**40).__ceil__())+1,10**20)
