"""Check exported files, numerical consistency and coordinate/component mapping."""
from pathlib import Path
import csv
import json
import numpy as np
from PIL import Image
from .pipeline import read_frame


def validate(root):
    out = Path(root)/'PZT_Atomiser_Analysis'
    manifest = json.loads((out/'results/figure_manifest.json').read_text())
    for stem in manifest:
        for suffix in ['.pdf','.svg','.png','.npz']:
            p = (out/stem).with_suffix(suffix)
            assert p.exists() and p.stat().st_size>0, p
        with Image.open((out/stem).with_suffix('.png')) as im:
            assert abs(im.info['dpi'][0]-600)<1
        data = np.load((out/stem).with_suffix('.npz'))
        if 'values' in data:
            z = data['values']
            lo,hi = data['color_limits']
            assert np.nanmin(z)>=lo-1e-9 and np.nanmax(z)<=hi+1e-9
    z = np.load(out/'results/Batch4_Cam1_statistics.npz')
    table = np.genfromtxt(out/'results/Batch4_Cam1_spatial_statistics.csv',delimiter=',',names=True)
    for key in table.dtype.names:
        np.testing.assert_allclose(table[key],z[key].ravel(),equal_nan=True)
    # Synthetic orthogonal reflection and scaling; paired mask must follow both components.
    cfg = dict(components=['u','v'],transform=[[1,0],[0,-1]],velocity_scale=.202,
               valid_flag_values=[1],flag_definition_source='synthetic test only')
    g = {'u':np.array([[2.,np.nan,5.]]),'v':np.array([[3.,4.,6.]]),'flag':np.array([[1,1,0]])}
    u,v = read_frame(g,cfg)
    np.testing.assert_allclose(u,[[.404,np.nan,np.nan]],equal_nan=True)
    np.testing.assert_allclose(v,[[-.606,np.nan,np.nan]],equal_nan=True)
    assert z['N'].min() == z['N'].max() == 1757
    report = {'figure_count':len(manifest),'formats':['pdf','svg','600-dpi png','npz'],
              'all_files_nonempty':True,'color_ranges_cover_data':True,
              'csv_matches_npz_all_rows':True,'synthetic_transform_and_pair_mask':True}
    (out/'results/validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    root=Path(__file__).resolve().parents[2]
    cfg=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text())
    if cfg.get('analysis_status') in ['provisional_calibrated','confirmed_calibrated']:
        from .validate_calibrated import main
        main()
    else:
        validate(root)
