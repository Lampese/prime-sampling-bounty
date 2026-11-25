from fractions import Fraction as F
import json

from . import ledger as v

ACCEPTED_JOINS = (F(651, 4), F(150))

def read_number(value):
    if isinstance(value,bool) or not isinstance(value,(int,str)):
        raise ValueError('Write endpoints and nonintegral numbers as strings or integers.')
    return F(value)

def grade(manifest, *, join=ACCEPTED_JOINS[0]):
    if join not in ACCEPTED_JOINS:
        raise ValueError('Choose an accepted coverage endpoint.')
    spec=json.loads(manifest.read_text())
    target=read_number(spec['target'])
    if not 40<=target<=join:
        raise ValueError(f'Choose a target in [40,{join}].')
    cells=[]
    start=target
    for cell in spec['cells']:
        left,right=map(read_number,(cell['left'],cell['right']))
        if left!=start or not left<right<=175:
            raise ValueError(f'Invalid coverage after {start}, next interval [{left},{right}].')
        cells.append((left,right,cell))
        start=right
    if start<join:
        raise ValueError(f'Coverage stops at {start}, before {join}.')
    minimum=None
    base=manifest.resolve().parent
    for left,right,cell in cells:
        tables=[]
        for item in cell.get('tables',[]):
            kind=item['kind']
            if kind not in ('single','pair','six'):
                raise ValueError(f'Unknown rate class {kind}.')
            tables.append(('extra',kind,left,right,int(item['grid']),base/item['file']))
        budgets=[base/path for path in cell.get('budgets',[])]
        lower=read_number(cell['main_lower']) if 'main_lower' in cell else None
        results=v.probe(left,right,tables,budgets,lower)
        for branch,result in results.items():
            if result[1]<=0:
                v.report(left,right,branch,result)
                raise ValueError(f'Strict positivity failed for {branch} on [{left},{right}].')
            minimum=result[1] if minimum is None else min(minimum,result[1])
    print(f'Continuous numerical coverage from {target} to {start}; archived coverage joins at {join}.')
    if minimum is not None:
        print(f'{len(cells)} submitted cells, both branches, minimum {v.fmt(minimum)}.')
