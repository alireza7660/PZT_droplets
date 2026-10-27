"""Validate final signed arrays/display masks without regenerating figures."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import nbformat as nb
from .final_analysis import validate,RENAME


def main():
    root=Path(__file__).resolve().parents[2]
    cfg=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text());out=root/cfg['output']
    with np.load(out/'data/statistics.npz') as z:
        fields={k:z[k].copy() for k in RENAME};r={'N':z['N'].copy(),'xy':np.stack([z['x_mm'],z['y_mm']])}
    manifest=json.loads((out/'methods/figure_manifest.json').read_text())
    validate(r,fields,cfg,out,SimpleNamespace(manifest=manifest))
    for entry in manifest:
        svg=(out/entry['files']['svg']).read_text(encoding='utf-8')
        assert not any(s in svg.lower() for s in ['provisional','uncalibrated','speed'])
    b=nb.read(out/'methods/PZT_final_analysis.ipynb',as_version=4)
    for cell in b.cells:
        if cell.cell_type=='code':
            assert cell.execution_count is not None
            assert not any(o.output_type=='error' for o in cell.outputs)
    print((out/'methods/summary.json').read_text())
    print('Final arrays, masks, exports and executed notebook verified.')


if __name__=='__main__':main()
