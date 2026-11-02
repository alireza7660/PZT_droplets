"""Verify the physical jet-under-plume entrance without a profile-crossing criterion."""
import json
from pathlib import Path
import numpy as np
import nbformat as nb
from PIL import Image,ImageDraw
from .aligned_codirected import load_context,prepare,pzt_sample
from .combined_flow_prediction import PROJECT,sha,write_json
from .droplet_dynamics import acceleration,reynolds,relaxation_time,drag_factor
from .prescribed_gas import UnsupportedPosition


def read(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))


def main():
    c=load_context();field,states,table,alignment,audit=prepare(c,write=False);out=c['out'];cfg=c['cfg'];props=c['props'];d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,props));Q=c['Q']
    np.testing.assert_allclose(Q@Q.T,np.eye(2),atol=1e-15);np.testing.assert_allclose(field.to_common(field.to_local([0,0])),[0,0],atol=1e-15);np.testing.assert_allclose(Q.T@np.array(cfg['gravity_vector_m_s2']),[0,-props['g_m_s2']],atol=1e-14)
    base=c['base'];p=c['profile'];valid=p['finite_pair']&np.isfinite(p['mean_u_p_m_s'])&np.isfinite(p['mean_v_p_m_s']);assert base['y_p_mm']==p['y_p_mm'][valid].min();assert not audit['crossing_used_to_select_entry']
    e=alignment['entry'];np.testing.assert_allclose(e['measured_Uc_m_s'],cfg['co_directed']['entrance_centreline_velocity_m_s'],rtol=0,atol=1e-12)
    bracket=alignment['developing_branch_bracket'];i,j=bracket['lower_native_column'],bracket['upper_native_column'];assert j==i+1 and c['valid'][i:j+1].all();np.testing.assert_allclose(bracket['x_j0_mm'],c['sx'][i]*(1-bracket['upper_weight'])+c['sx'][j]*bracket['upper_weight'],atol=1e-12)
    np.testing.assert_allclose(field.to_local([0,0]),[e['ZNMF_downstream_mm']*1e-3,e['ZNMF_local_centreline_Y_mm']*1e-3],atol=1e-15)
    gravity=np.array(cfg['gravity_vector_m_s2'])*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3']);ref=next(r for r in table if r['is_reference']);assert ref['droplet_v_p_m_s']>ref['gas_v_p_m_s'];assert ref['initial_drag_y_p_m_s2']<0
    np.testing.assert_allclose([ref['droplet_u_p_m_s'],ref['droplet_v_p_m_s']],[p['mean_u_p_m_s'][base['profile_index']],p['mean_v_p_m_s'][base['profile_index']]],atol=1e-12)
    saved=read(out/'data/initial_conditions.csv');assert len(saved)==len(states)
    for st,row in zip(states,table):
        initial=saved[saved['release_number']==st['number']][0]
        for key in ('h_mm','x_p_mm','droplet_u_p_m_s','droplet_v_p_m_s','gas_u_p_m_s','gas_v_p_m_s'):np.testing.assert_allclose(initial[key],row[key],atol=1e-12)
        assert initial['h_mm']==0
        q=read(out/'data'/f'release_{st["number"]:02d}.csv');points=np.c_[q['x_p_mm'],q['h_mm']]*1e-3;gas=np.array([field.sample(pos) for pos in points]);drop=np.c_[q['droplet_u_p_m_s'],q['droplet_v_p_m_s']];net=acceleration(drop,d,gas,props,cfg);drag=net-gravity;mass=q['droplet_mass_kg']
        np.testing.assert_allclose(gas,np.c_[q['gas_u_p_m_s'],q['gas_v_p_m_s']],atol=1e-13);np.testing.assert_allclose(drop[0],st['initial'][2:],atol=1e-15);np.testing.assert_allclose(q['Re_p'],reynolds(d,drop,gas,props),rtol=1e-12,atol=1e-13)
        np.testing.assert_allclose(q['ZNMF_downstream_mm'],e['ZNMF_downstream_mm']+q['h_mm'],atol=1e-12);np.testing.assert_allclose(q['y_p_mm'],q['h_mm']+base['y_p_mm'],atol=1e-12)
        lab=np.c_[q['x_p_mm'],q['h_mm']+base['y_p_mm']]@Q+c['origin'];np.testing.assert_allclose(lab,np.c_[q['PZT_lab_frame_x_mm'],q['PZT_lab_frame_y_mm']],atol=1e-12)
        np.testing.assert_allclose(drag,(drag_factor(q['Re_p'],cfg)/tau)[:,None]*(gas-drop),rtol=1e-10,atol=1e-9);assert np.all(np.sum(drag*(gas-drop),axis=1)>=-1e-12)
        for k,component in enumerate(('x_p','y_p')):
            np.testing.assert_allclose(q['net_a_'+component+'_m_s2'],net[:,k],rtol=1e-10,atol=1e-9);np.testing.assert_allclose(q['effective_gravity_a_'+component+'_m_s2'],gravity[k],atol=1e-14)
            residual=mass*(drop[:,k]-drop[0,k])-q['drag_impulse_'+component+'_N_s']-q['effective_gravity_impulse_'+component+'_N_s'];assert np.max(abs(residual/mass))<1e-10
            assert np.max(abs(q['force_balance_residual_'+component+'_m_s2']))<cfg['momentum_audit']['force_balance_tolerance_m_s2'];assert np.max(abs(q['momentum_residual_trapezoid_'+component+'_N_s']/mass))<cfg['momentum_audit']['independent_trapezoid_velocity_residual_tolerance_m_s']
        np.testing.assert_allclose(q['Stk_D_local'],tau*np.linalg.norm(gas,axis=1)/c['D0'],rtol=1e-12)
    for check in json.loads((out/'validation/convergence.json').read_text()):
        err=check['maximum_errors_x_y_u_v'];assert max(err[:2])<=cfg['solver']['position_convergence_tolerance_m'];assert max(err[2:])<=cfg['solver']['velocity_convergence_tolerance_m_s'];assert check['termination']==check['half_step_termination']
    candidates=read(out/'data/release_candidate_audit.csv');assert int(np.sum(candidates['accepted']))==len(states);assert len(candidates)==len(cfg['co_directed']['release_column_offsets'])
    for candidate in audit['release_candidates']:
        if candidate['accepted']:continue
        try:pzt_sample(c,[candidate['x_p_mm'],base['y_p_mm']]);field.sample([candidate['x_p_mm']*1e-3,0])
        except (ValueError,UnsupportedPosition):pass
        else:raise AssertionError('Rejected candidate is supported')
    try:field.sample(field.to_common([field.x.min()-1e-3,field.y.min()]))
    except UnsupportedPosition:pass
    else:raise AssertionError('Gas extrapolated outside measured support')
    events=json.loads((out/'data/reference_event_summary.json').read_text());minimum=events['first_velocity_minimum'];zero=events['first_axial_zero_slip']
    if minimum:assert abs(minimum['net_a_y_p_m_s2'])<1e-6 and minimum['slip_v_p_m_s']>0
    if zero:assert abs(zero['slip_v_p_m_s'])<1e-9;np.testing.assert_allclose(zero['net_a_y_p_m_s2'],gravity[1],atol=1e-6)
    result=json.loads((out/'results_summary.json').read_text());assert result['common_interval_h_mm'][0]==0
    force_tests=json.loads((out/'validation/analytical_force_checks.json').read_text());assert force_tests['zero_slip_drag_zero'] and force_tests['analytical_Stokes_constant_gas_2D_position_and_velocity_pass']
    exported=np.load(out/'data/prescribed_gas_plume_frame.npz');xx,hh,ug,vg=field.common_mesh();np.testing.assert_allclose(exported['h_mm'],hh*1e3,atol=1e-12);np.testing.assert_allclose(exported['mean_gas_v_p_m_s'],vg,atol=1e-12);np.testing.assert_array_equal(exported['finite_pair'],field.valid)
    protected=json.loads((out/'validation/authoritative_hashes.json').read_text());assert all(sha(p)==h for p,h in protected.items())
    manifest=json.loads((out/'report_manifest.json').read_text());paths=[]
    for stem in manifest['figures']:
        path=out/'figures'/stem
        for ext in ('.pdf','.svg','.png'):assert path.with_suffix(ext).stat().st_size>1000
        with Image.open(path.with_suffix('.png')) as im:assert min(im.info['dpi'])>599
        paths.append(path.with_suffix('.png'))
    for start in range(0,len(paths),6):
        canvas=Image.new('RGB',(1800,1100),'white');draw=ImageDraw.Draw(canvas)
        for k,path in enumerate(paths[start:start+6]):
            with Image.open(path) as im:im.thumbnail((590,510));x=(k%3)*600;y=(k//3)*550;canvas.paste(im,(x,y+24));draw.text((x+3,y+3),path.stem,fill='black')
        canvas.save(out/'validation'/f'visual_review_{start//6+1}.png')
    book=nb.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4);cells=[cell for cell in book.cells if cell.cell_type=='code' and cell.metadata.get('co_directed')];assert len(cells)==2 and all(cell.execution_count is not None for cell in cells);assert sum('image/png' in o.get('data',{}) for cell in cells for o in cell.outputs)==len(paths)
    receipt=dict(entrance_is_first_supported_measured_axis_point=True,centreline_target_interpolation_verified=True,crossing_not_used=True,initial_droplet_faster_than_sampled_gas=True,initial_axial_drag_negative=True,supported_candidate_count=len(states),rejected_candidate_count=len(candidates)-len(states),single_rigid_placement_verified=True,lab_plume_entrance_inverse_maps_verified=True,gravity_and_SI_units_verified=True,initial_native_components_and_weights_verified=True,slip_direction_force_and_momentum_verified=True,time_step_and_event_refinement_passed=True,unsupported_starts_rejected_without_extrapolation=True,all_saved_gas_samples_supported=True,authoritative_files_hash_verified=len(protected),Plume_Aligned_files_unchanged=sum('Plume_Aligned' in Path(p).parts for p in protected),executed_notebook_verified=True,figure_count=len(paths),PDF_SVG_600dpi_PNG_verified=True)
    write_json(out/'validation/completion_verification.json',receipt);print(json.dumps(receipt,indent=2));return receipt


if __name__=='__main__':main()
