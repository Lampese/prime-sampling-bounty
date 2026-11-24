"""Emit a single strip table with rational rounding. One table per invocation."""
from argparse import ArgumentParser
from fractions import Fraction as F

from ..rectangles.generate import emit
from ..spectrum.modes import MODES


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('left', type=F)
    parser.add_argument('right', type=F)
    parser.add_argument('--mode', choices=MODES, default='pair10')
    parser.add_argument('--grid', type=int, default=3200)
    args = parser.parse_args()
    if not __debug__ or not (150 <= args.left < args.right <= 175):
        parser.error('Use assertions and an interval contained in [150,175].')
    if args.grid < 600 or args.grid % 200:
        parser.error('Choose a positive multiple of 200, at least 600.')
    emit(args.left, args.right, args.mode, args.grid)


if __name__ == '__main__':
    main()
