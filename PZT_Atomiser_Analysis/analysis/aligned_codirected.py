"""Co-directed developed-droplet response using preserved plume-aligned measurements."""
import copy,json,shutil
from pathlib import Path
from datetime import datetime
import numpy as np
from .combined_flow_prediction import ROOT,PROJECT,sources_and_arrays,write_csv,write_json,sha
from .interaction_plane import bracket,supported_value
from .prescribed_gas import MeasuredGasField,integrate_droplet,trajectory_quantities,UnsupportedPosition
from .droplet_dynamics import relaxation_time,acceleration,reynolds,drag_factor


def read_npz(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}


def archive_existing(out):
    if not out.exists():return None
    dest=out.parent/'Archive'/('Co_Directed_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    assert dest.resolve().is_relative_to(PROJECT.resolve()) and out.resolve().is_relative_to(PROJECT.resolve())
    dest.mkdir(parents=True);shutil.move(str(out),str(dest/'generated'))
    receipt=dict(reason='Superseded generated co-directed prediction; measurements unchanged',generated_package=str((dest/'generated').relative_to(ROOT)))
    write_json(dest/'archive_manifest.json',receipt)
    return str((dest/'archive_manifest.json').relative_to(ROOT))


def load_context():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text())
    if cfg.get('interaction_work_status')=='paused_for_plume_alignment_review':raise RuntimeError('Interaction work is paused')
    laboratory,jet,xj,yj,D0,sources=sources_and_arrays(cfg)
    aligned_dir=ROOT/cfg['plume_aligned_directory'];a=read_npz(aligned_dir/'data/statistics_plume_aligned.npz');profile=read_npz(aligned_dir/'data/fitted_axis_axial_profile.npz')
    meta=json.loads((aligned_dir/'transformation.json').read_text());Q=np.asarray(meta['orthonormal_matrix_lab_to_plume']);origin=np.asarray(meta['reference_point_lab_mm'])
    # No calibration/reprocessing: confirm saved transforms against authoritative SI fields.
    positions=np.stack([laboratory['x_mm'],laboratory['y_mm']],-1);vectors=np.stack([laboratory['mean_u'],laboratory['mean_v']],-1)
    np.testing.assert_allclose((positions-origin)@Q.T,np.stack([a['x_p_mm'],a['y_p_mm']],-1),atol=1e-13)
    np.testing.assert_allclose(vectors@Q.T,np.stack([a['mean_u_p'],a['mean_v_p']],-1),atol=1e-14)
    np.testing.assert_array_equal(a['N'],laboratory['N'])
    np.testing.assert_allclose(cfg['gravity_vector_m_s2'],meta['gravity_plume_m_s2'],atol=1e-14)
    np.testing.assert_allclose(Q.T@np.array(cfg['gravity_vector_m_s2']),meta['gravity_lab_m_s2'],atol=1e-14)
    parent=ROOT/cfg['output'];out=parent/'Co_Directed';props=json.loads((parent/'fluid_properties.json').read_text())
    gas_cfg=json.loads((ROOT/cfg['znmf_configuration']).read_text());gas_path=(ROOT/cfg['znmf_configuration']).resolve().parent.parent/'tables/FOV_1_2_3_centreline_and_width.csv'
    gas=np.genfromtxt(gas_path,names=True,delimiter=',',dtype=None,encoding='utf-8')
    sx=gas['x_over_D0']*float(gas_cfg['diameter_mm']);uc=gas['Uc_m_per_s'];valid=(gas['status']=='measured')&np.isfinite(uc)&np.isfinite(gas['yc_absolute_over_D0'])
    for role,path in [('plume_aligned_native_statistics',aligned_dir/'data/statistics_plume_aligned.npz'),('fitted_axis_profile',aligned_dir/'data/fitted_axis_axial_profile.npz'),('plume_transformation',aligned_dir/'transformation.json'),('ZNMF_centreline',gas_path),('saved_fluid_properties',parent/'fluid_properties.json')]:sources['sources'].append(dict(role=role,path=str(path.resolve()),sha256=sha(path)))
    sources['pzt'].update(coordinates='Saved full-support x_p,y_p; laboratory coordinates and inverse Q retained',velocities='Both saved aligned mean_u_p,mean_v_p, already m/s; no calibration reapplied',profile='Fixed fitted axis x_p=0; old fixed laboratory-column profile excluded from model')
    sources['placement']=dict(common_frame=cfg['common_frame'],pzt_transform=meta['formulas'],gravity_plume_m_s2=cfg['gravity_vector_m_s2'],warning='Jet-under-plume design placement at first measured plane; profile crossing is not used')
    valid_axis=np.flatnonzero(profile['finite_pair']&np.isfinite(profile['mean_u_p_m_s'])&np.isfinite(profile['mean_v_p_m_s']))
    index=int(valid_axis[np.argmin(profile['y_p_mm'][valid_axis])]);base=dict(profile_index=index,y_p_mm=float(profile['y_p_mm'][index]),lab_x_mm=float(profile['lab_x_mm'][index]),lab_y_mm=float(profile['lab_y_mm'][index]),reference='First supported measured fitted-axis point, not documented physical atomiser surface')
    profile['h_mm']=profile['y_p_mm']-base['y_p_mm']
    return dict(base=base,cfg=cfg,lab=laboratory,a=a,profile=profile,meta=meta,Q=Q,origin=origin,jet=jet,xj=xj,yj=yj,D0=D0,sources=sources,props=props,out=out,gas=gas,gas_cfg=gas_cfg,gas_path=gas_path,sx=sx,uc=uc,valid=valid)


def pzt_sample(ctx,point_plume_mm):
    """Weights on original native rectangular cells, applied to actual aligned components."""
    point=np.asarray(point_plume_mm);lab=point@ctx['Q']+ctx['origin'];z=ctx['lab'];a=ctx['a']
    x=z['x_mm'][0];order=np.argsort(z['y_mm'][:,0]);y=z['y_mm'][order,0]
    j,wx=bracket(x,lab[0],np.ones(len(x),bool));i,wy=bracket(y,lab[1],np.ones(len(y),bool));rows=order[i:i+2]
    w=np.array([[(1-wx)*(1-wy),wx*(1-wy)],[(1-wx)*wy,wx*wy]])
    u=a['mean_u_p'][rows,j:j+2];v=a['mean_v_p'][rows,j:j+2]
    if not (np.isfinite(u).all() and np.isfinite(v).all()):raise ValueError('PZT interpolation cell contains a missing paired mean')
    uv=np.array([np.sum(w*u),np.sum(w*v)])
    labu=z['mean_u'][rows,j:j+2];labv=z['mean_v'][rows,j:j+2]
    np.testing.assert_allclose(uv,np.array([np.sum(w*labu),np.sum(w*labv)])@ctx['Q'].T,atol=1e-14)
    np.testing.assert_allclose(w.sum(),1.,atol=1e-15)
    return uv,dict(point_plume_mm=point.tolist(),point_lab_mm=lab.tolist(),native_rows_in_weight_order=rows.tolist(),native_columns=[j,j+1],weights=w.tolist(),source_mean_u_p=u.tolist(),source_mean_v_p=v.tolist(),source_mean_lab_u=labu.tolist(),source_mean_lab_v=labv.tolist(),source_N=z['N'][rows,j:j+2].tolist())


def prepare(ctx,write=True):
    from .base_placement import prepare_base
    return prepare_base(ctx,write)


def run():
    ctx=load_context();out=ctx['out']
    previous_summary=json.loads((out/'results_summary.json').read_text()) if (out/'results_summary.json').exists() else None
    archive=archive_existing(out)
    for name in ('data','figures','validation'):(out/name).mkdir(parents=True,exist_ok=True)
    cfg=copy.deepcopy(ctx['cfg']);cfg['solver']['maximum_time_s']=cfg['co_directed']['maximum_time_s'];cfg['solver']['time_step_s']=cfg['co_directed']['time_step_s'];ctx['cfg']=cfg
    protected={s['path']:s['sha256'] for s in ctx['sources']['sources']}
    final=out.parent.parent
    for path in final.rglob('*'):
        if path.is_file() and out.parent not in path.resolve().parents and ('Plume_Aligned' in path.parts or path.suffix.lower() in ('.npz','.csv','.pdf','.png','.svg')):protected[str(path.resolve())]=sha(path)
    field,states,initial,alignment,audit=prepare(ctx)
    from .momentum_audit import release_coverage,export_placed_gas
    export_placed_gas(ctx,field)
    coverage=release_coverage(ctx,states);write_json(out/'data/release_coverage.json',coverage)
    d=float(cfg['droplet_diameters_um'][0])*1e-6;assert len(cfg['droplet_diameters_um'])==1
    tau=float(relaxation_time(d,ctx['props']));histories=[];records=[];convergence=[]
    gravity=np.array(cfg['gravity_vector_m_s2'])*(1-ctx['props']['rho_air_kg_m3']/ctx['props']['rho_water_kg_m3'])
    for state in states:
        controls=copy.deepcopy(cfg)
        for refinement in range(6):
            tr,reason=integrate_droplet(field,state['initial'],d,ctx['props'],controls);fine,fine_reason=integrate_droplet(field,state['initial'],d,ctx['props'],controls,.5)
            if len(tr)<2 or len(fine)<2:raise ValueError('No supported trajectory after release')
            end=min(tr[-1,0],fine[-1,0]);times=np.linspace(0,end,401)
            errors=[float(np.max(abs(np.interp(times,tr[:,0],tr[:,k])-np.interp(times,fine[:,0],fine[:,k])))) for k in range(1,5)]
            if max(errors[:2])<=cfg['solver']['position_convergence_tolerance_m'] and max(errors[2:])<=cfg['solver']['velocity_convergence_tolerance_m_s']:break
            controls['solver']['time_step_s']*=.5
        else:raise AssertionError(('Step convergence failed',errors))
        q=trajectory_quantities(field,tr,d,ctx['props'],ctx['D0'],cfg);q['x_p_mm']=q.pop('x_mm');q['h_mm']=q.pop('y_mm');q['y_p_mm']=q['h_mm']+ctx['base']['y_p_mm']
        for old,new in [('droplet_u_m_s','droplet_u_p_m_s'),('droplet_v_m_s','droplet_v_p_m_s'),('gas_u_m_s','gas_u_p_m_s'),('gas_v_m_s','gas_v_p_m_s'),('slip_u_m_s','slip_u_p_m_s'),('slip_v_m_s','slip_v_p_m_s'),('Re_in_plane','Re_p'),('St_D_local','Stk_D_local')]:q[new]=q.pop(old)
        acc=acceleration(tr[:,3:5],d,np.c_[q['gas_u_p_m_s'],q['gas_v_p_m_s']],ctx['props'],cfg)
        q['net_a_x_p_m_s2']=acc[:,0];q['net_a_y_p_m_s2']=acc[:,1];q['drag_a_y_p_m_s2']=acc[:,1]-gravity[1]
        from .momentum_audit import add_force_and_impulse_fields, budget_record
        add_force_and_impulse_fields(q,tr,d,ctx['props'],cfg,acc)
        q['ZNMF_downstream_mm']=q['h_mm']-alignment['entry']['h_mm']+alignment['entry']['ZNMF_downstream_mm']
        uc_on=[]
        for s in q['ZNMF_downstream_mm']:
            try:uc_on.append(supported_value(ctx['sx'],ctx['uc'],ctx['valid'],s))
            except ValueError:uc_on.append(np.nan)
        q['measured_centreline_Uc_m_s']=np.array(uc_on)
        q['delta_x_p_mm']=q['x_p_mm']-q['x_p_mm'][0];q['delta_h_mm']=q['h_mm']-q['h_mm'][0]
        lab=np.c_[q['x_p_mm'],q['h_mm']+ctx['base']['y_p_mm']]@ctx['Q']+ctx['origin'];q['PZT_lab_frame_x_mm']=lab[:,0];q['PZT_lab_frame_y_mm']=lab[:,1]
        positive=np.flatnonzero(q['droplet_v_p_m_s']<=0);ascent_stop=int(positive[0])+1 if len(positive) else len(tr);peak=int(np.argmax(q['h_mm'][:ascent_stop]))
        number=state['number'];write_csv(out/'data'/f'release_{number:02d}.csv',q)
        h=dict(number=number,reference=state['reference'],q=q,peak=peak,termination=reason,fine_trajectory=fine if state['reference'] else None);histories.append(h)
        records.append(dict(release_number=number,is_reference=state['reference'],initial_x_p_mm=float(q['x_p_mm'][0]),initial_h_mm=float(q['h_mm'][0]),initial_u_p_m_s=float(q['droplet_u_p_m_s'][0]),initial_v_p_m_s=float(q['droplet_v_p_m_s'][0]),initial_gas_u_p_m_s=float(q['gas_u_p_m_s'][0]),initial_gas_v_p_m_s=float(q['gas_v_p_m_s'][0]),end_time_s=float(q['t_s'][-1]),end_h_mm=float(q['h_mm'][-1]),end_v_p_m_s=float(q['droplet_v_p_m_s'][-1]),delta_x_p_mm=float(q['delta_x_p_mm'][-1]),delta_v_p_m_s=float(q['droplet_v_p_m_s'][-1]-q['droplet_v_p_m_s'][0]),max_Re_p=float(np.max(q['Re_p'])),termination=reason))
        convergence.append(dict(release_number=number,time_step_s=controls['solver']['time_step_s'],maximum_errors_x_y_u_v=errors,termination=reason,half_step_termination=fine_reason,terminal_time_difference_s=float(abs(tr[-1,0]-fine[-1,0])),maximum_RK4_impulse_closure_equivalent_m_s=float(np.max(abs((tr[:,3:5]-state['initial'][2:4])-tr[:,5:7]-tr[:,7:9]))),maximum_half_step_RK4_impulse_closure_equivalent_m_s=float(np.max(abs((fine[:,3:5]-state['initial'][2:4])-fine[:,5:7]-fine[:,7:9])))))
    lower=alignment['entry']['h_mm'];upper=min(float(ctx['profile']['h_mm'][np.isfinite(ctx['profile']['mean_v_p_m_s'])].max()),min(float(h['q']['h_mm'][h['peak']]) for h in histories))
    if upper<=lower:raise ValueError('No common ascending interval')
    common=[]
    for h in histories:
        q=h['q'];n=h['peak']+1;yy=q['h_mm'][:n]
        vals={k:float(np.interp(upper,yy,q[k][:n])) for k in ('t_s','droplet_v_p_m_s','gas_v_p_m_s','x_p_mm')}
        baseline=supported_value(ctx['profile']['h_mm'],ctx['profile']['mean_v_p_m_s'],ctx['profile']['finite_pair'],upper)
        common.append(dict(release_number=h['number'],is_reference=h['reference'],lower_h_mm=lower,upper_h_mm=upper,transit_time_ms=vals['t_s']*1e3,entry_droplet_v_p_m_s=float(q['droplet_v_p_m_s'][0]),end_droplet_v_p_m_s=vals['droplet_v_p_m_s'],delta_v_p_m_s=vals['droplet_v_p_m_s']-float(q['droplet_v_p_m_s'][0]),end_sampled_gas_v_p_m_s=vals['gas_v_p_m_s'],PZT_only_axis_v_at_upper_m_s=baseline,lateral_displacement_mm=vals['x_p_mm']-float(q['x_p_mm'][0])))
    for name,rows in [('trajectory_summary',records),('common_height_comparison',common)]:write_csv(out/'data'/(name+'.csv'),{k:[r[k] for r in rows] for k in rows[0]})
    from .momentum_audit import write_momentum_tables,synthetic_force_checks
    momentum_summary=write_momentum_tables(ctx,histories,lower,upper)
    synthetic=synthetic_force_checks(cfg,ctx['props'],d);write_json(out/'validation/analytical_force_checks.json',synthetic)
    from .base_prediction_report import make_figures_and_report
    result=dict(reference_full_support=next(r for r in records if r['is_reference']),momentum=momentum_summary,release_coverage=coverage,successfully_integrated_cases=len(histories),placement_mode='plume_base_target_centreline',entry=alignment['entry'],common_interval_h_mm=[lower,upper],common_comparison=common,reference_comparison=next(r for r in common if r['is_reference']),tau_stokes_s=tau,nominal_diameter_um=d*1e6,gravity_plume_m_s2=cfg['gravity_vector_m_s2'],termination_reasons={r:sum(h['termination']==r for h in histories) for r in set(h['termination'] for h in histories)},archive_manifest=archive or cfg['co_directed']['archive_manifest'])
    if previous_summary:
        comparison={k:dict(previous=float(previous_summary['reference_comparison'][k]),current=float(result['reference_comparison'][k]),difference=float(result['reference_comparison'][k]-previous_summary['reference_comparison'][k])) for k in ('entry_droplet_v_p_m_s','end_droplet_v_p_m_s','transit_time_ms','lateral_displacement_mm')}
        write_json(out/'validation/previous_reference_comparison.json',comparison)
    from .trajectory_events import analyse_reference_events
    result['reference_events']=analyse_reference_events(ctx,field,histories,upper)
    stems=make_figures_and_report(ctx,field,histories,alignment,result)
    write_json(out/'results_summary.json',result);write_json(out/'configuration_used.json',cfg);write_json(out/'validation/convergence.json',convergence)
    code_paths=set((PROJECT/'analysis').glob('*codirected*.py'))|{PROJECT/'analysis'/name for name in ('co_directed.py','co_directed_report.py','droplet_dynamics.py','prescribed_gas.py','update_combined_flow_notebook.py','validate_co_directed.py','combined_flow_prediction.py','interaction_plane.py','prospective_studies.py','momentum_audit.py','base_placement.py','base_prediction_report.py','trajectory_events.py','validate_base_prediction.py')}
    ctx['sources']['reusable_code']={str(p):sha(p) for p in sorted(code_paths)}
    ctx['sources']['run_configuration']=dict(path=str(PROJECT/'config/combined_flow_prediction.json'),sha256=sha(PROJECT/'config/combined_flow_prediction.json'));ctx['sources']['aligned_placement']=alignment;ctx['sources']['properties']=ctx['props'];write_json(out/'source_manifest.json',ctx['sources'])
    write_json(out/'report_manifest.json',dict(figures=stems,numerical_files=[str(p.relative_to(out)) for p in (out/'data').glob('*')],text=['README.md','methods_results.tex','captions.tex'],entrance_target_interpolated_not_hardcoded_station=True))
    for path,digest in protected.items():assert sha(path)==digest,path
    write_json(out/'validation/authoritative_hashes.json',protected)
    write_json(out.parent/'prediction_status.json',dict(status='aligned_co_directed_prediction_complete',prediction_case_count=len(histories),active_package='Co_Directed',combined_experiment_conducted=False))
    write_json(out.parent/'report_manifest.json',dict(active_package='Co_Directed',manifest='Co_Directed/report_manifest.json',plume_frame='x_p,h; h=y_p-y_p_base',obsolete_prediction_displays_excluded=True))
    (out.parent/'README.md').write_text('Active prediction: Co_Directed/README.md. Jet-under-plume placement: the developing Uc=0.96 m/s station is at the first measured plume-axis plane h=0; profile crossings do not select entry. Uses full-support plume-aligned measured PZT fields. Reproduce from project root: .\\.venv\\Scripts\\python.exe -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook\n',encoding='utf-8')
    return result
