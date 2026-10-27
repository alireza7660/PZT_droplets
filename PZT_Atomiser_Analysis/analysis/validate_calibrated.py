"""Validate saved calibrated outputs without regenerating figures."""
import json
from pathlib import Path
import numpy as np
from .calibrated import LABELS, validate_outputs


def main():
    root=Path(__file__).resolve().parents[2]
    cfg=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text())
    out=root/cfg['output'];results={}
    for condition in cfg['conditions']:
        cid=condition['id']
        with np.load(out/'results'/f'{cid}_statistics.npz') as z:
            r=dict(fields={k:z[k].copy() for k in LABELS},N=z['N'].copy())
        with np.load(out/'results'/f'{cid}_station_samples.npz') as z:
            r['samples']=z['converted_V_U'].copy()
        results[cid]=r
    validate_outputs(out,results,cfg)
    print((out/'results/validation.json').read_text())


if __name__=='__main__':main()
