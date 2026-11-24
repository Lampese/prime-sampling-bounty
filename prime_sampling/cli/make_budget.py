"""Rebuild the four scalar budget entries on a closed interval."""
from argparse import ArgumentParser
from fractions import Fraction as F
from math import prod

from ..budget.build import build
from ..budget.certificate import eu, ps


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('left', type=F)
    p.add_argument('right', type=F)
    p.add_argument('--rank', type=int)
    p.add_argument('--method', choices=('multi', 'single'), default='multi')
    p.add_argument('--size', type=F)
    p.add_argument('--fine', type=F)
    args = p.parse_args()
    if not __debug__ or not 150 <= args.left < args.right <= 175:
        p.error('Use assertions and an interval inside [150,175].')
    rank = args.rank or next(i for i in range(1, 44) if prod(ps[:i+1]) > eu(args.right))
    if not 1 <= rank <= 43:
        p.error('Rank must lie in [1,43].')
    size = args.size or (F(3,40) if args.method == 'multi' else F(1,8))
    fine = args.fine or (F(1,1000) if args.method == 'multi' else F(1,500))
    if min(size, fine) <= 0:
        p.error('Positive mesh sizes are required.')
    print(args.left, args.right, rank, *build(args.left, args.right, rank, args.method, size, fine))


if __name__ == '__main__':
    main()
