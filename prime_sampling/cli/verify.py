"""Exact accounting of continuous parameter intervals. Python 3.10+."""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path

from ..accounting.frontier import verify_frontier
from ..accounting.ledger import probe, report


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('--detail', action='store_true')
    p.add_argument('--frontier', action='store_true', help='show both evaluations at the left endpoint')
    p.add_argument('--probe', type=F, nargs=2, metavar=('LEFT', 'RIGHT'))
    p.add_argument('--add', nargs=3, action='append', default=[], metavar=('KIND', 'GRID', 'FILE'),
                   help='extra rectangle table for the specified probe cell')
    p.add_argument('--budget', type=Path, action='append', default=[],
                   help='additional budget rows from the make_budget command for the probe')
    p.add_argument('--main-lower', type=F,
                   help='additional proven lower bound for A on the probe interval, in real units')
    args = p.parse_args()
    if not __debug__:
        p.error('Enable Python assertions.')
    if args.probe:
        a, z = args.probe
        if not 40 <= a < z <= 175:
            p.error('Choose a probe inside [40,175].')
        extras = []
        try:
            for kind, grid, file in args.add:
                if kind not in ('single', 'pair', 'six'):
                    p.error('KIND must be single, pair or six.')
                extras.append(('extra', kind, a, z, int(grid), Path(file)))
            evaluations = probe(a, z, extras, args.budget, args.main_lower)
        except ValueError as error:
            p.error(str(error))
        for branch, result in evaluations.items():
            report(a, z, branch, result)
        return
    if args.add or args.budget or args.main_lower is not None:
        p.error('Extra inputs require --probe.')
    verify_frontier(detail=args.detail, show_frontier=args.frontier)


if __name__ == '__main__':
    main()
