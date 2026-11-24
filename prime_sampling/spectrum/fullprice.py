from math import isqrt
from fractions import Fraction as F

from ..rectangles.lattice import dirichlet_cap,ceildiv
from ..rectangles.near import weighted_near
S=10**20
def root(x):return F(isqrt((x*S*S).__ceil__())+1,S)
WHEELS=((1,1),(2,1),(3,2),(6,2),(15,8),(30,8),(210,48),(2310,480))
