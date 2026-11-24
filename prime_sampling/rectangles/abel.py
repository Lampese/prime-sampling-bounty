from fractions import Fraction as F

def upper_mass(u,v,m):
    eta=m.thetaerr(min(u*m.zl,(1-v)*m.zl-m.ALOG))
    hh=(1+eta)*m.b(v/u,'ln',True)+2*eta/(u*m.zl)
    a=max(m.LAM,m.Xlo(2*u-1))
    A0=1-a+eta*(1+a)
    A1=(1+eta)*(1-a)-(1-eta)*a*m.b(1/a,'ln',False)
    assert A1>=0
    value=2*hh*(A0+A1/((1-v)*m.zl-m.ALOG))/(1-v)
    return F((value*10**22).__ceil__(),10**22)

def lower_count(cl,et,eps,z,m):
    assert 0<cl<et<=1
    low=(1-eps)*et-(1+eps)*cl
    correction=(1-eps)*(et-cl+et*(m.b(1/et,'ln',False) if et<1 else F(0)))-(1+eps)*cl*(m.b(1/cl,'ln',True) if cl<1 else F(0))
    value=low/z+max(F(0),correction)/z**2 if low>=0 else low/z
    return F((value*10**22).__floor__(),10**22)
