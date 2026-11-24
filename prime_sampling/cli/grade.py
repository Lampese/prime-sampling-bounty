"""Check continuous coverage of a submitted ledger, using both modulus branches."""
from argparse import ArgumentParser
from pathlib import Path

from ..accounting.submission import grade


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    args = p.parse_args()
    if not __debug__:
        p.error('Enable assertions.')
    try:
        grade(args.manifest)
    except (ValueError, KeyError, OSError, AssertionError) as error:
        p.exit(1, f'Ledger rejected: {error}\n')


if __name__ == '__main__':
    main()
