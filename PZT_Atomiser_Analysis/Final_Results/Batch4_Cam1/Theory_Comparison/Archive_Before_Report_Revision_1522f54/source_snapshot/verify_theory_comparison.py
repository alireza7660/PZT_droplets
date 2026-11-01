"""Validate saved model CSVs, exports, notebook and unchanged experiment files."""
import hashlib
import json
import numpy as np
import nbformat
from PIL import Image
from .theory_comparison import FINAL,PROJECT,reynolds,relaxation_time,drag_factor


def main():
    out=FINAL/'Theory_Comparison'
    record=json.loads((out/'parameters_used.json').read_text());p=record['properties'];cfg=record['configuration']
    hashes=json.loads((out/'experimental_input_hashes.json').read_text())
    for name,digest in hashes.items():assert hashlib.sha256((FINAL/name).read_bytes()).hexdigest()==digest,name
    measured=np.genfromtxt(out/'measured_fixed_x_profile.csv',names=True,delimiter=',')
    y0=record['initial_condition']['height_reference_y0_mm']
    np.testing.assert_allclose(measured['height_above_first_measurement_mm'],measured['y_mm']-y0,rtol=0,atol=1e-14)
    assert measured['height_above_first_measurement_mm'][0]==0
    with np.load(FINAL/'data/statistics.npz') as z:
        row=measured['native_row'].astype(int);col=measured['native_column'].astype(int)
        for key in ['mean_v','mean_u','N','y_mm','x_mm']:
            np.testing.assert_array_equal(measured[key],z[key][row,col])
    rc=np.genfromtxt(out/'reynolds_vs_height.csv',names=True,delimiter=',')
    dc=np.genfromtxt(out/'diameter_sensitivity.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    for i,d_um in enumerate(cfg['diameters_um']):
        d=d_um*1e-6;tau=relaxation_time(d,p)
        expected=reynolds(d,np.c_[measured['mean_u'],measured['mean_v']],cfg['gas_velocity_m_s'],p)
        np.testing.assert_allclose(rc[f'Re_still_air_d{d_um:g}um'],expected,rtol=1e-14)
        tr=np.genfromtxt(out/f'trajectory_quiescent_d{d_um:g}um.csv',names=True,delimiter=',')
        assert np.all(tr['v_m_s']>=0) and np.all(np.diff(tr['y_mm'])>=0)
        np.testing.assert_allclose(tr['height_above_first_measurement_mm'],tr['y_mm']-y0,rtol=0,atol=1e-14)
        assert abs(tr['height_above_first_measurement_mm'][0])<1e-14
        np.testing.assert_allclose(tr['tau_drag_secant_s'],tau/drag_factor(tr['Re_in_plane'],cfg),rtol=1e-14)
        np.testing.assert_allclose(dc['ascending_rise_within_domain_mm'][i],tr['y_mm'][-1]-tr['y_mm'][0],atol=1e-14)
    names=['velocity_theory_comparison','reynolds_vs_height','particle_relaxation_time','diameter_sensitivity','znmf_local_nozzle_stokes']
    for name in names:
        for extension in ['pdf','svg','png']:assert (out/(name+'.'+extension)).stat().st_size>1000
        with Image.open(out/(name+'.png')) as im:assert min(im.info['dpi'])>599
    book=nbformat.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4)
    cells=[c for c in book.cells if c.cell_type=='code' and c.metadata.get('theory_comparison')]
    assert len(cells)==1 and cells[0].execution_count is not None
    assert not any(o.output_type=='error' for o in cells[0].outputs)
    assert sum('image/png' in o.get('data',{}) for o in cells[0].outputs)==5
    result=dict(experimental_files_hash_verified=len(hashes),profile_CSV_exact_native_samples=True,relative_height_zero_and_alignment=True,
                Re_CSV_recalculated=True,trajectory_CSV_ascending_only=True,secant_time_CSV_verified=True,
                figure_count=5,all_three_formats=True,png_dpi=600,executed_notebook_figures=5,
                PDF_raster_review='Unavailable: no Poppler or PDFium renderer; exported PNGs visually reviewed',
                pdf_skill_marker='Node.js unavailable; marker invocation failed before authoring')
    (out/'completion_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
