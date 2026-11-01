"""Independent saved-output checks for the revised scientific theory package."""
import ast
import hashlib
import json
import numpy as np
import nbformat
from PIL import Image
from .theory_comparison import FINAL,PROJECT,reynolds,relaxation_time,drag_factor
from .theory_package import FIGURES,read_csv,digest


def main():
    out=FINAL/'Theory_Comparison';checkdir=out/'validation'
    params=json.loads((out/'parameters_used.json').read_text());p=params['properties'];cfg=params['configuration']
    state=params['initial_condition'];y0=state['initial_y_mm'];top=state['measurement_height_extent_mm']
    hashes=json.loads((checkdir/'experimental_input_hashes.json').read_text())
    for name,sha in hashes.items():assert digest(FINAL/name)==sha,name
    measured=read_csv(out/'measured_fixed_x_profile.csv')
    np.testing.assert_array_equal(measured['height_above_first_measurement_mm'],measured['y_mm']-y0)
    with np.load(FINAL/'data/statistics.npz') as z:
        row=measured['native_row'].astype(int);col=measured['native_column'].astype(int)
        for key in ['mean_v','mean_u','N','y_mm','x_mm']:np.testing.assert_array_equal(measured[key],z[key][row,col])
    re_mean=read_csv(out/'reynolds_measured_velocity_assumed_air.csv')
    ref_all=read_csv(out/'reynolds_reference_trajectories.csv')
    endpoints=np.genfromtxt(out/'reference_endpoints.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    old=out/'Archive_Before_Report_Revision_1522f54'
    unchanged_model_arrays=0
    for d_um in cfg['diameters_um']:
        d=d_um*1e-6;tau=relaxation_time(d,p)
        expected=reynolds(d,np.c_[measured['mean_u'],measured['mean_v']],cfg['gas_velocity_m_s'],p)
        np.testing.assert_array_equal(re_mean[f'Re_mean_vector_d{d_um:g}um'],expected)
        tr=read_csv(out/f'trajectory_stationary_air_d{d_um:g}um.csv')
        selection=ref_all[ref_all['assumed_diameter_um']==d_um]
        for key in tr.dtype.names:np.testing.assert_array_equal(selection[key],tr[key])
        assert np.all(tr['v_m_s']>=0) and np.all(np.diff(tr['height_above_first_measurement_mm'])>=0)
        np.testing.assert_allclose(tr['height_above_first_measurement_mm'],tr['y_mm']-y0,rtol=0,atol=1e-14)
        np.testing.assert_array_equal(tr['Re_in_plane'],reynolds(d,np.c_[tr['u_m_s'],tr['v_m_s']],cfg['gas_velocity_m_s'],p))
        np.testing.assert_allclose(tr['tau_drag_secant_s'],tau/drag_factor(tr['Re_in_plane'],cfg),rtol=1e-14)
        endpoint=endpoints[(endpoints['reference']=='stationary_air')&(endpoints['assumed_diameter_um']==d_um)][0]
        np.testing.assert_allclose(tr['height_above_first_measurement_mm'][-1],endpoint['plotted_endpoint_height_mm'],atol=1e-14)
        if endpoint['plotted_endpoint_reason']=='turning_point':
            assert tr['v_m_s'][-1]==0
            np.testing.assert_allclose(endpoint['actual_turning_height_mm'],endpoint['plotted_endpoint_height_mm'],atol=1e-14)
        elif endpoint['plotted_endpoint_reason']=='domain_top':
            np.testing.assert_allclose(tr['height_above_first_measurement_mm'][-1],top,atol=1e-12)
            assert endpoint['actual_turning_height_mm']>top
        else:raise AssertionError('Trajectory ended without a physical or domain endpoint')
        archived=read_csv(old/f'trajectory_quiescent_d{d_um:g}um.csv')
        for key in archived.dtype.names:np.testing.assert_array_equal(tr[key],archived[key])
        unchanged_model_arrays+=1
    gravity=endpoints[endpoints['reference']=='gravity_only'][0]
    np.testing.assert_allclose(gravity['actual_turning_height_mm'],state['initial_v_m_s']**2/(2*p['g_m_s2'])*1000,rtol=1e-14)
    np.testing.assert_allclose(gravity['plotted_endpoint_height_mm'],min(top,gravity['actual_turning_height_mm']),rtol=1e-14)
    ballistic=read_csv(out/'ballistic_reference.csv')
    np.testing.assert_allclose(ballistic['v_m_s']**2,state['initial_v_m_s']**2-2*p['g_m_s2']*ballistic['height_above_first_measurement_mm']*.001,atol=1e-14)
    old_ast=ast.parse((old/'source_snapshot/theory_comparison.py').read_text(encoding='utf-8'))
    new_ast=ast.parse((PROJECT/'analysis/theory_comparison.py').read_text(encoding='utf-8'))
    functions=['fluid_properties','relaxation_time','reynolds','drag_factor','acceleration','upward_trajectory','validate_physics']
    for name in functions:
        a=next(n for n in old_ast.body if isinstance(n,ast.FunctionDef) and n.name==name)
        b=next(n for n in new_ast.body if isinstance(n,ast.FunctionDef) and n.name==name)
        assert ast.dump(a)==ast.dump(b),name
    manifest=json.loads((out/'figure_curve_manifest.json').read_text())
    assert [m['stem'] for m in manifest]==list(FIGURES)
    curves_verified=0
    for item in manifest:
        for c in item['curves']:
            data=read_csv(out/c['csv'])
            if c['filter_column']:data=data[data[c['filter_column']]==c['filter_value']]
            for dimension in ('x','y'):
                arr=np.asarray(data[c[dimension+'_column']]*c[dimension+'_factor'],dtype='<f8')
                assert hashlib.sha256(arr.tobytes()).hexdigest()==c[dimension+'_sha256']
            assert len(data)==c['point_count'];curves_verified+=1
        for ext in ('pdf','svg','png'):assert (out/(item['stem']+'.'+ext)).stat().st_size>1000
        with Image.open(out/(item['stem']+'.png')) as im:assert min(im.info['dpi'])>599
    assert manifest[0]['curves']==manifest[1]['curves']
    assert manifest[1]['xlim']==[0.,2.]
    assert manifest[0]['ylim']==manifest[1]['ylim']
    assert manifest[-1]['yscale']=='log'
    jet=read_csv(out/'Jet_Interaction_Sensitivity/znmf_local_nozzle_stokes.csv')
    for d_um in cfg['diameters_um']:
        np.testing.assert_array_equal(jet[f'St_D_d{d_um:g}um'],relaxation_time(d_um*1e-6,p)*jet['Uc_m_s']/cfg['znmf_nozzle_diameter_m'])
    tau_table=read_csv(out/'particle_relaxation_time.csv')
    np.testing.assert_array_equal(tau_table['tau_stokes_s'],relaxation_time(tau_table['assumed_diameter_um']*1e-6,p))
    report=json.loads((out/'report_manifest.json').read_text())
    assert all('Archive' not in f['pdf'] and (out/f['pdf']).exists() for f in report['figures'])
    for obsolete in ['velocity_theory_comparison.pdf','reynolds_vs_height.csv','znmf_local_nozzle_stokes.pdf']:assert not (out/obsolete).exists()
    book=nbformat.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4)
    cells=[c for c in book.cells if c.cell_type=='code' and c.metadata.get('theory_comparison')]
    assert len(cells)==1 and cells[0].execution_count is not None
    assert not any(o.output_type=='error' for o in cells[0].outputs)
    assert sum('image/png' in o.get('data',{}) for o in cells[0].outputs)==len(FIGURES)
    result=dict(experimental_files_hash_verified=len(hashes),profile_exact_native_samples=True,
        full_near_velocity_curve_arrays_and_styles_identical=True,mean_and_trajectory_Re_separate=True,
        actual_turning_and_truncated_endpoints_separate=True,physics_source_functions_unchanged=len(functions),
        original_trajectory_arrays_unchanged=unchanged_model_arrays,figure_CSV_curves_verified=curves_verified,
        jet_Stokes_definition_verified=True,jet_log_ordinate=True,archive_excluded_from_report=True,
        figure_count=len(FIGURES),png_dpi=600,executed_notebook_figures=len(FIGURES),
        PDF_raster_review='No PDF renderer available; final exported PNGs visually reviewed')
    (checkdir/'saved_output_checks.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
