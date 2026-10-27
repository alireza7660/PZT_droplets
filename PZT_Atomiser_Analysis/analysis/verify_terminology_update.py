"""Verify terminology-only numerical preservation and remove superseded files."""
from pathlib import Path
import json
import numpy as np
import nbformat as nb
from PIL import Image


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'provisional_calibrated'
    baseline=out/'diagnostics/terminology_baseline'
    checks=[]
    for prior in baseline.glob('*_statistics.npz'):
        with np.load(prior) as old,np.load(out/'results'/prior.name) as new:
            for key in old.files:
                newkey=key.replace('speed','velocity_magnitude')
                np.testing.assert_array_equal(old[key],new[newkey])
                checks.append(newkey)
    manifest=json.loads((out/'results/figure_manifest.json').read_text())
    previous=json.loads((baseline/'prior_figure_manifest.json').read_text())
    removed=[]
    for stem in previous:
        replacement=stem.replace('speed','velocity_magnitude')
        if replacement==stem:
            continue
        assert replacement in manifest
        for suffix in ['.pdf','.svg','.png','.npz']:
            p=(out/stem).with_suffix(suffix).resolve()
            q=(out/replacement).with_suffix(suffix).resolve()
            assert p.is_relative_to(out.resolve()) and q.exists()
            if p.exists():p.unlink();removed.append(str(p.relative_to(out)))
    for stem in manifest:
        assert 'speed' not in stem
        svg=(out/stem).with_suffix('.svg').read_text(encoding='utf-8')
        assert 'speed' not in svg.lower()
    for path in [out/'METHODS.md',out/'CAPTIONS.md',out/'CAPTIONS.tex',out/'FINDINGS.md',root/'README.md']:
        assert 'speed' not in path.read_text(encoding='utf-8').lower(),path
    b=nb.read(root/'notebooks/Spatial_maps.ipynb',as_version=4)
    for cell in b.cells:
        if cell.cell_type=='code':
            assert cell.execution_count is not None
            assert not any(o.output_type=='error' for o in cell.outputs)
        assert 'speed' not in cell.source.lower()
    previews=out/'diagnostics/visual_review';previews.mkdir(exist_ok=True)
    for name in ['mean_U','profile_mean_U','U_pdf','frame_879_quiver']:
        with Image.open(out/'figures/Batch4_Cam1'/(name+'.png')) as im:
            im.thumbnail((950,800));im.save(previews/(name+'.png'))
    report=dict(numerical_arrays_bitwise_unchanged=checks,figure_count=len(manifest),
                superseded_files_removed=removed,figure_and_report_terminology_checked=True,
                notebook_executed_without_errors=True,
                coordinate_labels='Axial coordinate x; transverse coordinate y')
    (out/'results/terminology_update_verification.json').write_text(json.dumps(report,indent=2))
    print('Verified numerical preservation, figure labels, documentation and notebook; removed',len(removed),'superseded files.')


if __name__=='__main__':main()
