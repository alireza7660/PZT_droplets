"""Synthetic solver tests and preservation checks; no experimental size/placement cases."""
import copy
import json
import numpy as np
import nbformat
from PIL import Image
from .combined_flow_prediction import ROOT,PROJECT,write_json,sha,missing_inputs,sources_and_arrays,release_state
from .droplet_dynamics import fluid_properties,acceleration,relaxation_time,reynolds,drag_factor
from .prescribed_gas import MeasuredGasField,UnsupportedPosition,integrate_droplet,trajectory_quantities


def synthetic_checks(cfg):
    """Manufactured fields/states are test fixtures, never physical predictions."""
    p=fluid_properties(cfg);x=np.linspace(-1,1,9);y=np.linspace(-1,1,9);xx,yy=np.meshgrid(x,y)
    u=1+2*xx+3*yy;v=-.2+xx-yy
    field=MeasuredGasField(x,y,u,v,[.4,-.7],90)
    point=np.array([.12,-.13]);physical=field.to_common(point)
    np.testing.assert_allclose(field.to_local(physical),point,atol=1e-15)
    np.testing.assert_allclose(field.sample(physical),field.vectors_to_common([1+2*point[0]+3*point[1],-.2+point[0]-point[1]]),atol=1e-14)
    # Rotation acts identically on gas position and velocity; gravity remains common-frame downward.
    np.testing.assert_allclose(field.vectors_to_common([1,0]),[0,1],atol=1e-15)
    d=5e-6;zero=np.zeros(2);g=p['g_m_s2']*(1-p['rho_air_kg_m3']/p['rho_water_kg_m3'])
    np.testing.assert_allclose(acceleration(zero,d,zero,p,cfg),[0,-g])
    assert acceleration([0,1],d,zero,p,cfg)[1]<-g
    assert acceleration([0,-1],d,zero,p,cfg)[1]>-g
    assert acceleration(zero,d,[1,0],p,cfg)[0]>0
    np.testing.assert_allclose(reynolds(d,[3,4],zero,p),p['rho_air_kg_m3']*d*5/p['mu_air_Pa_s'])
    tiny=np.array([1e-9,2e-9]);a=acceleration(tiny,d,zero,p,cfg)-[0,-g]
    np.testing.assert_allclose(a,-tiny/relaxation_time(d,p),rtol=1e-6,atol=1e-14)
    assert float(drag_factor(0,cfg))==1
    try:drag_factor(800,cfg)
    except ValueError:pass
    else:raise AssertionError('Invalid drag Reynolds range accepted')
    quiet=MeasuredGasField(x,y,np.zeros_like(xx),np.zeros_like(xx),[0,0],0)
    controls=copy.deepcopy(cfg);controls['solver']['maximum_time_s']=.0002
    initial=np.array([0,0,.001,0.]);tau=float(relaxation_time(d,p))
    tr,reason=integrate_droplet(quiet,initial,d,p,controls);t=tr[:,0];decay=np.exp(-t/tau)
    np.testing.assert_allclose(tr[:,3],initial[2]*decay,rtol=1e-6,atol=1e-12)
    np.testing.assert_allclose(tr[:,4],-g*tau*(1-decay),rtol=1e-6,atol=1e-12)
    np.testing.assert_allclose(tr[:,1],initial[2]*tau*(1-decay),rtol=1e-6,atol=1e-12)
    assert reason=='maximum_time' and tr[-1,4]<0  # no old upward-only stop
    fine,_=integrate_droplet(quiet,initial,d,p,controls,.5)
    np.testing.assert_allclose(tr[-1,1:],fine[-1,1:],rtol=1e-6,atol=1e-12)
    q=trajectory_quantities(quiet,tr,d,p,.002,cfg)
    assert np.all(q['St_D_local']==0) and np.isinf(q['tau_flow_local_D0_over_gas_magnitude_s']).all()
    # Source support and holes are not extrapolated or bridged.
    try:quiet.sample([2,0])
    except UnsupportedPosition:pass
    else:raise AssertionError('Gas extrapolation occurred')
    bad=np.zeros_like(xx);bad[4,4]=np.nan
    masked=MeasuredGasField(x,y,bad,np.zeros_like(xx),[0,0],0)
    try:masked.sample([.1,.1])
    except UnsupportedPosition:pass
    else:raise AssertionError('Interpolation crossed missing gas vector')
    bx=np.linspace(0,.001,9);by=np.linspace(-.001,.001,9);bu=np.ones((9,9));bv=np.zeros((9,9))
    bounded=MeasuredGasField(bx,by,bu,bv,[0,0],0)
    controls=copy.deepcopy(cfg);controls['solver']['maximum_time_s']=.002
    exit_tr,reason=integrate_droplet(bounded,[.0005,0,1,0],20e-6,p,controls)
    np.testing.assert_array_equal(exit_tr[:,3],np.ones(len(exit_tr)))
    assert reason=='outside_measured_jet_domain' and exit_tr[-1,1]<=bx[-1]
    assert bx[-1]-exit_tr[-1,1]<1e-8
    for pos in exit_tr[:,1:3]:bounded.sample(pos)
    return dict(synthetic_only=True,coordinate_and_vector_rotation=True,bilinear_affine_field_exact=True,
        gravity_downward=True,drag_opposes_vector_slip=True,zero_slip_finite=True,
        low_Re_analytic_Stokes_limit=True,half_step_convergence=True,finite_Re_validity_enforced=True,
        no_upward_only_stop=True,source_boundary_stopping=True,missing_gas_cell_rejected=True,
        zero_gas_local_Stokes_handled=True)


def main():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    checks=synthetic_checks(cfg);write_json(out/'validation/synthetic_solver_tests.json',checks)
    protected=json.loads((out/'validation/authoritative_file_hashes.json').read_text())
    for path,expected in protected.items():assert sha(path)==expected,path
    archive=json.loads((ROOT/cfg['archive_manifest']).read_text())
    for item in archive['records']:assert sha(ROOT/item['archive_path'])==item['sha256'],item['archive_path']
    pzt,gas,xj,yj,D0,source=sources_and_arrays(cfg)
    gas_csv=np.genfromtxt(out/'data/znmf_only_local_mean_field.csv',names=True,delimiter=',')
    np.testing.assert_array_equal(gas_csv['U_m_s'],gas['U'].ravel())
    np.testing.assert_array_equal(gas_csv['V_m_s'],gas['V'].ravel())
    # Audit the dimensionalisation without changing authoritative velocity components.
    with np.load(ROOT/cfg['znmf_field']) as z:
        np.testing.assert_array_equal(gas['U'],z['U']);np.testing.assert_array_equal(gas['V'],z['V'])
        np.testing.assert_array_equal(xj,z['x']*D0)
        np.testing.assert_array_equal(yj,(z['y']-source['znmf']['nozzle_y_registered_D0'])*D0)
    for component in ('u','v'):
        with np.load(out/f'data/pzt_only_mean_{component}.npz') as z:
            np.testing.assert_array_equal(z['values'],pzt['mean_'+component])
            np.testing.assert_array_equal(z['display_mask'],~np.isfinite(pzt['mean_'+component]))
    profile=np.genfromtxt(out/'data/pzt_only_fixed_line_profile.csv',names=True,delimiter=',')
    i,j=profile['native_row'].astype(int),profile['native_column'].astype(int)
    for key,field in [('mean_v','mean_v'),('mean_u','mean_u'),('y_mm','y_mm'),('x_mm','x_mm'),('N','N')]:
        np.testing.assert_array_equal(profile[key],pzt[field][i,j])
    report=json.loads((out/'report_manifest.json').read_text());status=json.loads((out/'prediction_status.json').read_text())
    if cfg['droplet_diameters_um']:
        diameter_table=np.atleast_1d(np.genfromtxt(out/'data/droplet_model_parameters.csv',names=True,delimiter=','))
        np.testing.assert_array_equal(diameter_table['configured_diameter_um'],cfg['droplet_diameters_um'])
        np.testing.assert_array_equal(diameter_table['diameter_m'],np.asarray(cfg['droplet_diameters_um'])*1e-6)
        np.testing.assert_array_equal(diameter_table['tau_stokes_s'],relaxation_time(diameter_table['diameter_m'],fluid_properties(cfg)))
    for item in report['figures']:
        assert 'Archive' not in item['stem'] and 'Theory_Comparison' not in item['stem']
        for path in item['exports']:assert (out/path).stat().st_size>1000
        with Image.open(out/'figures'/(item['stem']+'.png')) as im:assert min(im.info['dpi'])>599
    if missing_inputs(cfg):
        assert status['status']=='awaiting_prediction_inputs' and status['prediction_case_count']==0
        assert not list((out/'figures').glob('predicted_*'))
    book=nbformat.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4)
    assert not any(c.metadata.get('theory_comparison') for c in book.cells)
    cell=[c for c in book.cells if c.cell_type=='code' and c.metadata.get('combined_flow_prediction')]
    assert len(cell)==1 and not any(o.output_type=='error' for o in cell[0].outputs)
    assert sum('image/png' in o.get('data',{}) for o in cell[0].outputs)==len(report['figures'])
    receipt=dict(authoritative_files_hash_verified=len(protected),archived_files_hash_verified=len(archive['records']),
        measured_PZT_arrays_masks_profiles_unchanged=True,ZNMF_velocities_not_recalibrated=True,
        coordinate_dimensionalisation_verified=True,synthetic_tests_passed=True,diameters_explicitly_sourced=bool(cfg.get('diameter_source_or_assumption')),
        prediction_status=status['status'],figure_count=len(report['figures']),PNG_dpi=600,notebook_executed=True)
    write_json(out/'validation/completion_verification.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
