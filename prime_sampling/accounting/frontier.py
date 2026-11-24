from fractions import Fraction as F

from .catalog import PANELS
from .ledger import evaluate, fmt, report

STAGES = {
    'A': (F(169), F(175), 24),
    'B': (F(657, 4), F(169), 54),
    'C': (F(655, 4), F(657, 4), 16),
    'D': (F(327, 2), F(655, 4), 4),
    'E': (F(651, 4), F(327, 2), 20),
}


def verify_frontier(*, detail=False, show_frontier=False):
    panels = sorted((F(a), F(z), stage, ref) for a, z, stage, ref in PANELS)
    assert len(panels) == 118 and panels[0][0] == F(651, 4) and panels[-1][1] == 175
    assert all(a < z for a, z, *_ in panels)
    assert all(x[1] == y[0] for x, y in zip(panels, panels[1:]))
    for stage, (left, right, count) in STAGES.items():
        intervals = [row for row in panels if row[2] == stage]
        assert len(intervals) == count and (intervals[0][0], intervals[-1][1]) == (left, right)
        worst = None
        for a, z, _, ref in intervals:
            branches = ('with3', 'without3') if stage in ('C', 'E') else ('all',)
            assert set(ref) == set(branches)
            for branch in branches:
                result = evaluate(a, z, stage, branch)
                assert tuple(ref[branch]) == result[:3] and result[1] > 0, (a, z, branch, result, ref)
                worst = result[1] if worst is None else min(worst, result[1])
                if detail or (show_frontier and a == F(651, 4)):
                    report(a, z, branch, result)
        print(f'{stage}  [{left}, {right}]  {count} intervals  minimum {fmt(worst)}', flush=True)
    print('118 adjoining closed intervals, 154 branch evaluations; lower endpoint 651/4.')
