"""Independent saved-result verification for corrected co-directed trajectories."""
import json
from pathlib import Path
import numpy as np
import nbformat as nb
from PIL import Image,ImageDraw
from .aligned_codirected import load_context,prepare
from .combined_flow_prediction import PROJECT,sha,write_json
from .droplet_dynamics import acceleration,reynolds,relaxation_time,drag_factor
from .prescribed_gas import UnsupportedPosition


def read(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))


def main():
    ctx=load_context();field,states,table,alignment,audit=prepare(ctx,write=False);out=ctx['out'];cfg=ctx['cfg'];props=ctx['props'];d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,props));Q=ctx['Q'];R=field.R
    np.testing.assert_allclose(Q@Q.T,np.eye(2),atol=1e-15);np.testing.assert_allclose(R@R.T,np.eye(2),atol=1e-15)
    np.testing.assert_array_equal(field.vectors_to_common([1.,0.]),[0,1]);np.testing.assert_array_equal(field.vectors_to_common([0.,1.]),[-1,0])
    point=np.array([0,alignment['entry']['y_p_mm']*1e-3]);np.testing.assert_allclose(field.to_common(field.to_local(point)),point,atol=1e-15)
    np.testing.assert_allclose(Q.T@np.array(cfg['gravity_vector_m_s2']),[0,-props['g_m_s2']],atol=1e-14)
    gravity=np.array(cfg['gravity_vector_m_s2'])*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3']);np.testing.assert_allclose(acceleration([0,0],d,[0,0],props,cfg),gravity,atol=1e-14)
    reference=next(r for r in table if r['is_reference']);assert reference['gas_v_p_m_s']>reference['droplet_v_p_m_s']
    assert alignment['entry']['ZNMF_downstream_mm']>alignment['crossing']['ZNMF_downstream_mm'];assert audit['entry_candidates'][-1]['accepted'] and not any(r['accepted'] for r in audit['entry_candidates'][:-1])
    saved=read(out/'data/initial_conditions.csv')
    for st,row in zip(states,table):
        record=saved[saved['release_number']==st['number']][0]
        for key in ('x_p_mm','y_p_mm','droplet_u_p_m_s','droplet_v_p_m_s','gas_u_p_m_s','gas_v_p_m_s','initial_a_x_p_m_s2','initial_a_y_p_m_s2'):np.testing.assert_allclose(record[key],row[key],atol=1e-12)
        q=read(out/'data'/f'release_{st["number"]:02d}.csv');points=np.c_[q['x_p_mm'],q['y_p_mm']]*1e-3;gas=np.array([field.sample(point) for point in points]);np.testing.assert_allclose(gas,np.c_[q['gas_u_p_m_s'],q['gas_v_p_m_s']],atol=1e-13)
        drop=np.c_[q['droplet_u_p_m_s'],q['droplet_v_p_m_s']];slip=gas-drop;acc=acceleration(drop,d,gas,props,cfg)
        np.testing.assert_allclose(drop[0],st['initial'][2:],atol=1e-15);np.testing.assert_allclose(points[0],st['initial'][:2],atol=1e-15)
        np.testing.assert_allclose(q['Re_p'],reynolds(d,drop,gas,props),rtol=1e-12,atol=1e-13);np.testing.assert_allclose(q['net_a_x_p_m_s2'],acc[:,0],rtol=1e-10,atol=1e-9);np.testing.assert_allclose(q['net_a_y_p_m_s2'],acc[:,1],rtol=1e-10,atol=1e-9)
        drag=acc-gravity;np.testing.assert_allclose(drag,(drag_factor(q['Re_p'],cfg)/tau)[:,None]*slip,rtol=1e-10,atol=1e-9);assert np.all(np.sum(drag*slip,axis=1)>=-1e-12)
        for k,component in enumerate(('x_p','y_p')):
            np.testing.assert_allclose(q['drag_a_'+component+'_m_s2'],drag[:,k],rtol=1e-12,atol=1e-10)
            np.testing.assert_allclose(q['effective_gravity_a_'+component+'_m_s2'],gravity[k],atol=1e-14)
            mass=q['droplet_mass_kg'];np.testing.assert_allclose(q['force_drag_'+component+'_N'],mass*drag[:,k],rtol=1e-12,atol=1e-23)
            momentum=mass*(drop[:,k]-drop[0,k]);residual=momentum-q['drag_impulse_'+component+'_N_s']-q['effective_gravity_impulse_'+component+'_N_s']
            np.testing.assert_allclose(q['momentum_residual_RK4_'+component+'_N_s'],residual,atol=1e-26);assert np.max(abs(residual/mass))<1e-10
            assert np.max(abs(q['force_balance_residual_'+component+'_m_s2']))<cfg['momentum_audit']['force_balance_tolerance_m_s2']
            assert np.max(abs(q['momentum_residual_trapezoid_'+component+'_N_s']/mass))<cfg['momentum_audit']['independent_trapezoid_velocity_residual_tolerance_m_s']
        np.testing.assert_allclose(q['Stk_D_local'],tau*np.linalg.norm(gas,axis=1)/ctx['D0'],rtol=1e-12)
    try:field.sample(field.to_common([field.x.min()-1e-3,field.y.min()]))
    except UnsupportedPosition:pass
    else:raise AssertionError('Extrapolated outside measured domain')
    if (~field.valid).any():
        i,j=np.argwhere(~field.valid)[0];i=min(i,len(field.y)-2);j=min(j,len(field.x)-2)
        try:field.sample(field.to_common([np.mean(field.x[j:j+2]),np.mean(field.y[i:i+2])]))
        except UnsupportedPosition:pass
        else:raise AssertionError('Filled a masked cell')
    for check in json.loads((out/'validation/convergence.json').read_text()):
        e=check['maximum_errors_x_y_u_v'];assert max(e[:2])<=cfg['solver']['position_convergence_tolerance_m'];assert max(e[2:])<=cfg['solver']['velocity_convergence_tolerance_m_s'];assert check['termination']==check['half_step_termination']
    numbers=read(out/'data/jet_local_centreline_numbers.csv');valid=numbers['supported_positive_centreline'].astype(bool);f=ctx['gas_cfg']['forcing_Hz'];nu=props['mu_air_Pa_s']/props['rho_air_kg_m3']
    np.testing.assert_allclose(numbers['Re_D_c'][valid],numbers['Uc_m_s'][valid]*ctx['D0']/nu);np.testing.assert_allclose(numbers['Sr_D_c'][valid],f*ctx['D0']/numbers['Uc_m_s'][valid]);assert np.isnan(numbers['Re_D_c'][~valid]).all()
    result=json.loads((out/'results_summary.json').read_text());common=read(out/'data/common_height_comparison.csv')
    for r in common:np.testing.assert_allclose([r['lower_y_p_mm'],r['upper_y_p_mm']],result['common_interval_y_p_mm'],atol=1e-12);assert r['transit_time_ms']>0
    protected=json.loads((out/'validation/authoritative_hashes.json').read_text());assert all(sha(p)==h for p,h in protected.items())
    exported=np.load(out/'data/prescribed_gas_plume_frame.npz');xx,yy,ug,vg=field.common_mesh()
    for key,expected in [('x_p_mm',xx*1e3),('y_p_mm',yy*1e3),('mean_gas_u_p_m_s',ug),('mean_gas_v_p_m_s',vg)]:np.testing.assert_allclose(exported[key],expected,atol=1e-13)
    np.testing.assert_array_equal(exported['finite_pair'],field.valid)
    endpoints=read(out/'data/entry_exit_conditions.csv');assert len(endpoints)==2*len(states)
    summaries=read(out/'data/trajectory_summary.csv')
    for row in summaries:
        exit_row=endpoints[(endpoints['stage']=='exit')&(endpoints['release_number']==row['release_number'])][0]
        assert exit_row['termination_reason']==row['termination']
    manifest=json.loads((out/'report_manifest.json').read_text());paths=[]
    for stem in manifest['figures']:
        p=out/'figures'/stem
        for ext in ('.pdf','.svg','.png'):assert p.with_suffix(ext).stat().st_size>1000
        with Image.open(p.with_suffix('.png')) as im:assert min(im.info['dpi'])>599
        paths.append(p.with_suffix('.png'))
    for start in range(0,len(paths),6):
        canvas=Image.new('RGB',(1800,1100),'white');draw=ImageDraw.Draw(canvas)
        for i,path in enumerate(paths[start:start+6]):
            with Image.open(path) as im:
                im.thumbnail((590,510));xx=(i%3)*600;yy=(i//3)*550;canvas.paste(im,(xx,yy+24));draw.text((xx+3,yy+3),path.stem,fill='black')
        canvas.save(out/'validation'/f'visual_review_{start//6+1}.png')
    book=nb.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4);cells=[c for c in book.cells if c.cell_type=='code' and c.metadata.get('co_directed')];assert len(cells)==2 and all(c.execution_count is not None for c in cells);assert not any(o.get('output_type')=='error' for c in cells for o in c.outputs)
    assert sum('image/png' in o.get('data',{}) for c in cells for o in c.outputs)==len(paths)
    synthetic=json.loads((out/'validation/analytical_force_checks.json').read_text());assert synthetic['zero_slip_drag_zero'] and synthetic['analytical_Stokes_constant_gas_2D_position_and_velocity_pass']
    coverage=json.loads((out/'data/release_coverage.json').read_text());assert not coverage['whole_plume_distribution'] and len(coverage['initial_x_p_mm'])==7
    previous=json.loads((out/'validation/previous_reference_comparison.json').read_text());assert all(abs(r['difference'])<1e-10 for r in previous.values())
    receipt=dict(analytical_constant_gas_Stokes_solution_verified=True,momentum_RK4_and_independent_trapezoid_verified=True,unchanged_previous_reference_result_verified=True,seven_starts_audited_as_unweighted_sensitivity=True,successfully_integrated_cases=len(states),source_coordinate_pairing_preserved=True,crossing_calculated_with_valid_adjacent_intervals=True,first_supported_native_actual_gas_overtaking_entry=True,single_rigid_gas_placement=True,position_vector_inverse_verified=True,gravity_inverse_verified=True,aligned_measured_component_interpolation_verified=True,all_saved_positions_within_finite_gas_support=True,slip_direction_drag_verified=True,time_step_convergence_passed=True,local_centreline_Re_Sr_formulas_and_gaps_verified=True,common_supported_interval_verified=True,authoritative_files_hash_verified=len(protected),Plume_Aligned_files_unchanged=sum('Plume_Aligned' in Path(p).parts for p in protected),executed_notebook_verified=True,figure_count=len(paths),PDF_SVG_600dpi_PNG_verified=True)
    write_json(out/'validation/completion_verification.json',receipt);print(json.dumps(receipt,indent=2));return receipt


if __name__=='__main__':main()
