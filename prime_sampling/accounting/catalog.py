from fractions import Fraction as F
import json
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / 'data'
INDEX = json.loads((DATA_DIR / 'index.json').read_text(encoding='utf-8'))

PANELS = INDEX['panels']
TABLES = [
    (group, kind, F(left), F(right), grid, DATA_DIR / name)
    for group, kind, left, right, grid, name in INDEX['tables']
]
BUDGETS = [
    (group, F(left), F(right), rank, values)
    for group, left, right, rank, values in INDEX['budgets']
]
