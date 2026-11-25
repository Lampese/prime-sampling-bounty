"""Check continuous coverage of a submitted ledger, using both modulus branches."""
from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path

from ..accounting.submission import ACCEPTED_JOINS, grade


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('--join', type=Fraction, choices=ACCEPTED_JOINS, default=ACCEPTED_JOINS[0],
                   help='join an accepted proof at this endpoint (default: %(default)s)')
    args = p.parse_args()
    if not __debug__:
        p.error('Enable assertions.')
    try:
        grade(args.manifest, join=args.join)
    except (ValueError, KeyError, OSError, AssertionError) as error:
        p.exit(1, f'Ledger rejected: {error}\n')


if __name__ == '__main__':
    main()
