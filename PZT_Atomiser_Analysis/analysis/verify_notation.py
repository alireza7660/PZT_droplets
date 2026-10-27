"""Snapshot/check numerical preservation and concise mathematical labels."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
import nbformat as nb

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'provisional_calibrated'
BASE=OUT/'diagnostics/notation_baseline'
KEYS=['components','coordinate_scale','velocity_scale','component_mapping',
      'coordinate_mapping','transform','origin','valid_flag_values']


def snapshot():
    BASE.mkdir(exist_ok=True)
    for path in (OUT/'results').glob('*.npz'):shutil.copy2(path,BASE/path.name)
    cfg=json.loads((ROOT/'config/analysis.json').read_text())
    (BASE/'measurement_config.json').write_text(json.dumps({k:cfg[k] for k in KEYS},indent=2))
    print('Saved numerical baseline before regeneration.')


def verify():
    from .calibrated import LABELS,AXES
    checked=[]
    for path in BASE.glob('*.npz'):
        with np.load(path) as old,np.load(OUT/'results'/path.name) as new:
            assert old.files==new.files,path
            for key in old.files:np.testing.assert_array_equal(old[key],new[key])
            checked.append(path.name)
    cfg=json.loads((ROOT/'config/analysis.json').read_text())
    assert {k:cfg[k] for k in KEYS}==json.loads((BASE/'measurement_config.json').read_text())
    assert AXES==dict(xlabel=r'$y$ ($\mathrm{mm}$)',ylabel=r'$x$ ($\mathrm{mm}$)')
    assert r'\overline{U}' in LABELS['mean_U'] and r'\overline{V}' in LABELS['mean_V']
    assert r"\overline{u'v'}" in LABELS['covariance']
    manifest=json.loads((OUT/'results/figure_manifest.json').read_text())
    for stem in manifest:
        svg=(OUT/stem).with_suffix('.svg').read_text(encoding='utf-8')
        assert 'speed' not in svg.lower()
        assert not any(s in svg for s in ['Transverse coordinate','Axial coordinate','Transverse position','Axial position'])
        for ext in ['.pdf','.svg','.png','.npz']:assert (OUT/stem).with_suffix(ext).stat().st_size>0
    selection=json.loads((OUT/'results/report_figure_selection.json').read_text())
    assert [Path(s).name for s in selection['primary_velocity_contours']]==['mean_U','mean_V']
    assert all('magnitude' in s for s in selection['supplementary_velocity_magnitude_figures'])
    b=nb.read(ROOT/'notebooks/Spatial_maps.ipynb',as_version=4)
    for cell in b.cells:
        if cell.cell_type=='code':
            assert cell.execution_count is not None
            assert not any(o.output_type=='error' for o in cell.outputs)
    report={'numeric_result_files_bitwise_equal':checked,'measurement_config_unchanged':True,
            'figure_count':len(manifest),'mathematical_labels_verified':True,
            'primary_velocity_contours':['mean_U','mean_V'],'magnitude_figures_supplementary':True,
            'notebook_executed_without_errors':True}
    (OUT/'results/notation_verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--snapshot',action='store_true')
    args=p.parse_args();snapshot() if args.snapshot else verify()
