"""Co-directed developed-droplet response using preserved plume-aligned measurements."""
import copy,json,shutil
from pathlib import Path
from datetime import datetime
import numpy as np
from .combined_flow_prediction import ROOT,PROJECT,sources_and_arrays,write_csv,write_json,sha
from .interaction_plane import bracket,supported_value,COMPARISON
from .prescribed_gas import MeasuredGasField,integrate_droplet,trajectory_quantities,UnsupportedPosition
from .droplet_dynamics import relaxation_time,acceleration,reynolds,drag_factor


def read_npz(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}


def archive_existing(out):
    if not out.exists():return None
    dest=out.parent/'Archive'/('Co_Directed_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    assert dest.resolve().is_relative_to(PROJECT.resolve()) and out.resolve().is_relative_to(PROJECT.resolve())
    dest.mkdir(parents=True);shutil.move(str(out),str(dest/'generated'))
    (dest/'source').mkdir()
    for path in (PROJECT/'analysis').glob('*codirected*.py'):shutil.copy2(path,dest/'source'/path.name)
    shutil.copy2(PROJECT/'config/combined_flow_prediction.json',dest/'source/configuration.json')
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
    for role,path in [('plume_aligned_native_statistics',aligned_dir/'data/statistics_plume_aligned.npz'),('fitted_axis_profile',aligned_dir/'data/fitted_axis_axial_profile.npz'),('plume_transformation',aligned_dir/'transformation.json'),('ZNMF_centreline',gas_path),('saved_fluid_properties',parent/'fluid_properties.json'),('original_comparison_script',COMPARISON/'plot_comparison.py'),('original_comparison_verification',COMPARISON/'verification.json')]:sources['sources'].append(dict(role=role,path=str(path.resolve()),sha256=sha(path)))
    sources['pzt'].update(coordinates='Saved full-support x_p,y_p; laboratory coordinates and inverse Q retained',velocities='Both saved aligned mean_u_p,mean_v_p, already m/s; no calibration reapplied',profile='Fixed fitted axis x_p=0; old fixed laboratory-column profile excluded from model')
    sources['placement']=dict(common_frame=cfg['common_frame'],pzt_transform=meta['formulas'],gravity_plume_m_s2=cfg['gravity_vector_m_s2'],warning='Original comparison coordinate pairing is a design placement, not proof of coincident physical origins')
    return dict(cfg=cfg,lab=laboratory,a=a,profile=profile,meta=meta,Q=Q,origin=origin,jet=jet,xj=xj,yj=yj,D0=D0,sources=sources,props=props,out=out,gas=gas,gas_cfg=gas_cfg,gas_path=gas_path,sx=sx,uc=uc,valid=valid)


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


def coordinates_at_pair(ctx,s):
    """Preserve the old comparison pairing in laboratory coordinates, not new axial zero."""
    height=s+ctx['cfg']['co_directed']['relative_axial_offset_mm'];fit=ctx['meta']['primary_fit'];lab=np.array([fit['intercept_x_at_lab_y0_mm']+fit['slope_dx_dy']*height,height]);p=(lab-ctx['origin'])@ctx['Q'].T
    return dict(ZNMF_downstream_mm=float(s),PZT_lab_x_mm=float(lab[0]),PZT_lab_y_mm=float(lab[1]),x_p_mm=float(p[0]),y_p_mm=float(p[1]))


def prepare(ctx,write=True):
    p=ctx['profile'];sx,uc,valid=ctx['sx'],ctx['uc'],ctx['valid'];py=p['lab_y_mm'];pv=p['mean_v_p_m_s'];pvalid=p['finite_pair']&np.isfinite(pv)&np.isfinite(p['mean_u_p_m_s']);offset=ctx['cfg']['co_directed']['relative_axial_offset_mm']
    # Original comparison plotted ZNMF x*D0 and original PZT laboratory y on 0..27 mm.
    original=json.loads((COMPARISON/'verification.json').read_text());domain=original['axis_limits_mm'];low=max(domain[0],float(sx[0]),float(py[0]-offset));high=min(domain[1],float(sx[-1]),float(py[-1]-offset))
    archive_gas=np.genfromtxt(COMPARISON/'znmf_plotted_data.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    archived_mask=(ctx['gas']['status']=='measured')&np.isfinite(uc)&(sx>=domain[0])&(sx<=domain[1])
    np.testing.assert_array_equal(sx[archived_mask],archive_gas['s_mm']);np.testing.assert_array_equal(uc[archived_mask],archive_gas['mean_axial_velocity_m_s'])
    nodes=np.unique(np.r_[low,high,sx[(sx>=low)&(sx<=high)],py[(py-offset>=low)&(py-offset<=high)]-offset]);roots=[]
    for left,right in zip(nodes[:-1],nodes[1:]):
        try:
            gi,_=bracket(sx,(left+right)/2,valid);pi,_=bracket(py,(left+right)/2+offset,pvalid)
            gv=np.interp([left,right],sx[gi:gi+2],uc[gi:gi+2]);dv=np.interp(np.array([left,right])+offset,py[pi:pi+2],pv[pi:pi+2])
        except ValueError:continue
        diff=gv-dv
        if diff[0]<=0<=diff[1] and diff[1]>diff[0]:
            root=float(left-diff[0]*(right-left)/(diff[1]-diff[0]))
            if roots and abs(roots[-1]['ZNMF_downstream_mm']-root)<1e-10:continue
            roots.append(dict(**coordinates_at_pair(ctx,root),common_axial_velocity_m_s=supported_value(sx,uc,valid,root),gas_native_lower_column=gi,PZT_profile_lower_index=pi,bracket_mm=[float(left),float(right)]))
    if not roots:raise ValueError('No supported gas-overtaking crossing')
    crossing=roots[0];attempts=[];selected=None
    nozzle_reference=json.loads((ROOT/ctx['cfg']['znmf_nozzle_reference']).read_text())['geometric_nozzle_centre_registered_D0']
    for j in np.flatnonzero((sx>crossing['ZNMF_downstream_mm'])&(sx<=high)&valid):
        coord=coordinates_at_pair(ctx,float(sx[j]));yc=(ctx['gas']['yc_absolute_over_D0'][j]-nozzle_reference)*ctx['D0']
        translation=np.array([yc,coord['y_p_mm']*1e-3-sx[j]*1e-3]);field=MeasuredGasField(ctx['xj'],ctx['yj'],ctx['jet']['U'],ctx['jet']['V'],translation,90.)
        try:
            vd,stencil=pzt_sample(ctx,[0,coord['y_p_mm']]);vg=field.sample(np.array([0,coord['y_p_mm']])*1e-3)
        except (ValueError,UnsupportedPosition) as problem:
            attempts.append(dict(ZNMF_native_column=int(j),**coord,accepted=False,reason=str(problem)));continue
        accepted=bool(vg[1]>vd[1]);attempts.append(dict(ZNMF_native_column=int(j),**coord,accepted=accepted,gas_v_p_m_s=float(vg[1]),PZT_v_p_m_s=float(vd[1]),slip_v_p_m_s=float(vg[1]-vd[1])))
        if accepted:selected=dict(**coord,ZNMF_native_column=int(j),ZNMF_local_centreline_Y_mm=float(yc*1e3),measured_Uc_m_s=float(uc[j]),sampled_gas_v_p_m_s=float(vg[1]),PZT_axis_v_p_m_s=float(vd[1]));break
    if selected is None:raise ValueError('No native downstream station has supported actual gas overtaking')
    np.testing.assert_allclose(field.to_local([0,selected['y_p_mm']*1e-3]),[selected['ZNMF_downstream_mm']*1e-3,selected['ZNMF_local_centreline_Y_mm']*1e-3],atol=1e-14)
    dx=float(np.median(np.diff(ctx['a']['x_p_mm'],axis=1)));states=[];table=[];stencils=[];d=ctx['cfg']['droplet_diameters_um'][0]*1e-6;props=ctx['props']
    gravity=np.array(ctx['cfg']['gravity_vector_m_s2'])*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3'])
    for offset_col in ctx['cfg']['co_directed']['release_column_offsets']:
        transverse=offset_col*dx;point=np.array([transverse,selected['y_p_mm']]);v,stencil=pzt_sample(ctx,point)
        try:gas=field.sample(point*1e-3)
        except UnsupportedPosition as problem:raise ValueError(f'Configured transverse release {transverse} mm unsupported: {problem}')
        state=np.r_[point*1e-3,v];acc=acceleration(v,d,gas,props,ctx['cfg']);drag=acc-gravity
        np.testing.assert_allclose(drag,(gas-v)*float(drag_factor(reynolds(d,v,gas,props),ctx['cfg'])/relaxation_time(d,props)),rtol=1e-12,atol=1e-12)
        assert np.dot(drag,gas-v)>=-1e-12
        local=field.to_local(point*1e-3);ix,wx=bracket(ctx['xj'],local[0],np.ones(len(ctx['xj']),bool));iy,wy=bracket(ctx['yj'],local[1],np.ones(len(ctx['yj']),bool))
        gw=np.array([[(1-wx)*(1-wy),wx*(1-wy)],[(1-wx)*wy,wx*wy]])
        record=dict(release_number=len(states)+1,is_reference=offset_col==0,transverse_grid_offset=int(offset_col),x_p_mm=float(point[0]),y_p_mm=float(point[1]),PZT_lab_x_mm=stencil['point_lab_mm'][0],PZT_lab_y_mm=stencil['point_lab_mm'][1],ZNMF_downstream_mm=float(local[0]*1e3),ZNMF_local_Y_mm=float(local[1]*1e3),PZT_native_row_0=stencil['native_rows_in_weight_order'][0],PZT_native_row_1=stencil['native_rows_in_weight_order'][1],PZT_native_col_0=stencil['native_columns'][0],PZT_native_col_1=stencil['native_columns'][1],PZT_weight_00=float(stencil['weights'][0][0]),PZT_weight_01=float(stencil['weights'][0][1]),PZT_weight_10=float(stencil['weights'][1][0]),PZT_weight_11=float(stencil['weights'][1][1]),droplet_u_p_m_s=float(v[0]),droplet_v_p_m_s=float(v[1]),gas_u_p_m_s=float(gas[0]),gas_v_p_m_s=float(gas[1]),measured_centreline_Uc_m_s=selected['measured_Uc_m_s'],sampled_minus_centreline_v_m_s=float(gas[1]-selected['measured_Uc_m_s']),slip_u_p_m_s=float(gas[0]-v[0]),slip_v_p_m_s=float(gas[1]-v[1]),initial_a_x_p_m_s2=float(acc[0]),initial_a_y_p_m_s2=float(acc[1]),initial_drag_x_p_m_s2=float(drag[0]),initial_drag_y_p_m_s2=float(drag[1]),Re_p=float(reynolds(d,v,gas,props)))
        stencil.update(release_number=len(states)+1,gas_lower_native_row=iy,gas_lower_native_column=ix,gas_weights=gw.tolist(),gas_U_corners=ctx['jet']['U'][iy:iy+2,ix:ix+2].tolist(),gas_V_corners=ctx['jet']['V'][iy:iy+2,ix:ix+2].tolist())
        states.append(dict(initial=state,reference=offset_col==0,number=len(states)+1));table.append(record);stencils.append(stencil)
    entry_yp=selected['y_p_mm'];entry_s=selected['ZNMF_downstream_mm']
    alignment=dict(entry=selected,crossing=crossing,gas_to_plume_rotation_matrix=field.R.tolist(),gas_nozzle_translation_plume_m=field.origin.tolist(),PZT_lab_to_plume_Q=ctx['Q'].tolist(),PZT_reference_lab_mm=ctx['origin'].tolist(),gravity_plume_m_s2=ctx['cfg']['gravity_vector_m_s2'],comparison_pairing='s_ZNMF=y_PZT_lab-offset; y_p=(y_PZT_lab-y_ref)/cos(theta) on fitted axis',rigid_placement='[x_p,y_p]=R90*[Xj,Yj]+translation; vectors [u_p,v_p]=[-Vj,Uj]. One placement only; no downstream recentering.',gas_streamwise_to_plume_formula=f'y_p_mm={entry_yp!r}+(ZNMF_downstream_mm-{entry_s!r})',interpretation='Prescribed design placement of independent datasets, not coincident experimental origins or physical atomiser relocation. Coordinate transform of measured inclined plume, not asserted physical remounting.')
    audit=dict(crossing=crossing,all_supported_overtaking_crossings=roots,selected_entry=selected,entry_candidates=attempts,original_comparison_domain_mm=domain,comparison_mapping_evidence=str(COMPARISON/'plot_comparison.py'),source_curve_changed='PZT fixed laboratory-column mean_v replaced by saved fixed fitted-axis mean_v_p; original laboratory-height pairing retained.',initial_conditions_exported_before_integration=True,calibration_reapplied=False,report_crop_used=False,gas_centreline_difference='Measured Uc is interpolated at each native measured centreline. A field sample follows the bilinear two-component field at the release/trajectory location. At the chosen native entry column the reference matches Uc to rounding; off-axis samples differ and downstream the fixed jet placement retains measured centreline wandering.')
    if write:
        out=ctx['out'];write_csv(out/'data/initial_conditions.csv',{k:[r[k] for r in table] for k in table[0]});write_json(out/'data/initial_interpolation_stencils.json',stencils);write_json(out/'interaction_plane_audit.json',audit);write_json(out/'transformation.json',alignment)
        write_csv(out/'data/PZT_corrected_axis_profile.csv',p)
        gy_pair=(sx+offset-ctx['origin'][1])/ctx['Q'][0,0];gy_rigid=entry_yp+sx-entry_s
        write_csv(out/'data/ZNMF_profile_coordinate_mapping.csv',dict(ZNMF_downstream_mm=sx,comparison_pair_PZT_lab_y_mm=sx+offset,comparison_y_p_mm=gy_pair,rigidly_placed_y_p_mm=gy_rigid,measured_Uc_m_s=uc,centreline_valid=valid,status=ctx['gas']['status']))
    return field,states,table,alignment,audit


def run():
    ctx=load_context();out=ctx['out'];archive=archive_existing(out)
    for name in ('data','figures','validation'):(out/name).mkdir(parents=True,exist_ok=True)
    cfg=copy.deepcopy(ctx['cfg']);cfg['solver']['maximum_time_s']=cfg['co_directed']['maximum_time_s'];cfg['solver']['time_step_s']=cfg['co_directed']['time_step_s'];ctx['cfg']=cfg
    protected={s['path']:s['sha256'] for s in ctx['sources']['sources']}
    final=out.parent.parent
    for path in final.rglob('*'):
        if path.is_file() and out.parent not in path.resolve().parents and ('Plume_Aligned' in path.parts or path.suffix.lower() in ('.npz','.csv','.pdf','.png','.svg')):protected[str(path.resolve())]=sha(path)
    field,states,initial,alignment,audit=prepare(ctx)
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
        q=trajectory_quantities(field,tr,d,ctx['props'],ctx['D0'],cfg);q['x_p_mm']=q.pop('x_mm');q['y_p_mm']=q.pop('y_mm')
        for old,new in [('droplet_u_m_s','droplet_u_p_m_s'),('droplet_v_m_s','droplet_v_p_m_s'),('gas_u_m_s','gas_u_p_m_s'),('gas_v_m_s','gas_v_p_m_s'),('slip_u_m_s','slip_u_p_m_s'),('slip_v_m_s','slip_v_p_m_s'),('Re_in_plane','Re_p'),('St_D_local','Stk_D_local')]:q[new]=q.pop(old)
        acc=acceleration(tr[:,3:5],d,np.c_[q['gas_u_p_m_s'],q['gas_v_p_m_s']],ctx['props'],cfg)
        q['net_a_x_p_m_s2']=acc[:,0];q['net_a_y_p_m_s2']=acc[:,1];q['drag_a_y_p_m_s2']=acc[:,1]-gravity[1]
        q['ZNMF_downstream_mm']=q['y_p_mm']-alignment['entry']['y_p_mm']+alignment['entry']['ZNMF_downstream_mm']
        uc_on=[]
        for s in q['ZNMF_downstream_mm']:
            try:uc_on.append(supported_value(ctx['sx'],ctx['uc'],ctx['valid'],s))
            except ValueError:uc_on.append(np.nan)
        q['measured_centreline_Uc_m_s']=np.array(uc_on)
        q['delta_x_p_mm']=q['x_p_mm']-q['x_p_mm'][0];q['delta_y_p_mm']=q['y_p_mm']-q['y_p_mm'][0]
        lab=np.c_[q['x_p_mm'],q['y_p_mm']]@ctx['Q']+ctx['origin'];q['PZT_lab_frame_x_mm']=lab[:,0];q['PZT_lab_frame_y_mm']=lab[:,1]
        positive=np.flatnonzero(q['droplet_v_p_m_s']<=0);ascent_stop=int(positive[0])+1 if len(positive) else len(tr);peak=int(np.argmax(q['y_p_mm'][:ascent_stop]))
        number=state['number'];write_csv(out/'data'/f'release_{number:02d}.csv',q)
        h=dict(number=number,reference=state['reference'],q=q,peak=peak,termination=reason);histories.append(h)
        records.append(dict(release_number=number,is_reference=state['reference'],initial_x_p_mm=float(q['x_p_mm'][0]),initial_y_p_mm=float(q['y_p_mm'][0]),initial_u_p_m_s=float(q['droplet_u_p_m_s'][0]),initial_v_p_m_s=float(q['droplet_v_p_m_s'][0]),initial_gas_u_p_m_s=float(q['gas_u_p_m_s'][0]),initial_gas_v_p_m_s=float(q['gas_v_p_m_s'][0]),end_time_s=float(q['t_s'][-1]),end_y_p_mm=float(q['y_p_mm'][-1]),end_v_p_m_s=float(q['droplet_v_p_m_s'][-1]),delta_x_p_mm=float(q['delta_x_p_mm'][-1]),delta_v_p_m_s=float(q['droplet_v_p_m_s'][-1]-q['droplet_v_p_m_s'][0]),max_Re_p=float(np.max(q['Re_p'])),termination=reason))
        convergence.append(dict(release_number=number,time_step_s=controls['solver']['time_step_s'],maximum_errors_x_y_u_v=errors,termination=reason,half_step_termination=fine_reason,terminal_time_difference_s=float(abs(tr[-1,0]-fine[-1,0]))))
    lower=alignment['entry']['y_p_mm'];upper=min(float(ctx['profile']['y_p_mm'][np.isfinite(ctx['profile']['mean_v_p_m_s'])].max()),min(float(h['q']['y_p_mm'][h['peak']]) for h in histories))
    if upper<=lower:raise ValueError('No common ascending interval')
    common=[]
    for h in histories:
        q=h['q'];n=h['peak']+1;yy=q['y_p_mm'][:n]
        vals={k:float(np.interp(upper,yy,q[k][:n])) for k in ('t_s','droplet_v_p_m_s','gas_v_p_m_s','x_p_mm')}
        baseline=supported_value(ctx['profile']['y_p_mm'],ctx['profile']['mean_v_p_m_s'],ctx['profile']['finite_pair'],upper)
        common.append(dict(release_number=h['number'],is_reference=h['reference'],lower_y_p_mm=lower,upper_y_p_mm=upper,transit_time_ms=vals['t_s']*1e3,entry_droplet_v_p_m_s=float(q['droplet_v_p_m_s'][0]),end_droplet_v_p_m_s=vals['droplet_v_p_m_s'],delta_v_p_m_s=vals['droplet_v_p_m_s']-float(q['droplet_v_p_m_s'][0]),end_sampled_gas_v_p_m_s=vals['gas_v_p_m_s'],PZT_only_axis_v_at_upper_m_s=baseline,lateral_displacement_mm=vals['x_p_mm']-float(q['x_p_mm'][0])))
    for name,rows in [('trajectory_summary',records),('common_height_comparison',common)]:write_csv(out/'data'/(name+'.csv'),{k:[r[k] for r in rows] for k in rows[0]})
    from .aligned_codirected_report import make_figures_and_report
    result=dict(successfully_integrated_cases=len(histories),crossing=alignment['crossing'],entry=alignment['entry'],common_interval_y_p_mm=[lower,upper],common_comparison=common,reference_comparison=next(r for r in common if r['is_reference']),tau_stokes_s=tau,nominal_diameter_um=d*1e6,gravity_plume_m_s2=cfg['gravity_vector_m_s2'],termination_reasons={r:sum(h['termination']==r for h in histories) for r in set(h['termination'] for h in histories)},archive_manifest=archive or cfg['co_directed']['archive_manifest'])
    stems=make_figures_and_report(ctx,field,histories,alignment,result)
    write_json(out/'results_summary.json',result);write_json(out/'configuration_used.json',cfg);write_json(out/'validation/convergence.json',convergence)
    ctx['sources']['reusable_code']={str(p):sha(p) for p in (PROJECT/'analysis').glob('*codirected*.py')};ctx['sources']['aligned_placement']=alignment;ctx['sources']['properties']=ctx['props'];write_json(out/'source_manifest.json',ctx['sources'])
    write_json(out/'report_manifest.json',dict(figures=stems,numerical_files=[str(p.relative_to(out)) for p in (out/'data').glob('*')],text=['README.md','methods_results.tex','captions.tex'],interaction_entry_calculated_not_hardcoded=True))
    for path,digest in protected.items():assert sha(path)==digest,path
    write_json(out/'validation/authoritative_hashes.json',protected)
    write_json(out.parent/'prediction_status.json',dict(status='aligned_co_directed_prediction_complete',prediction_case_count=len(histories),active_package='Co_Directed',combined_experiment_conducted=False))
    write_json(out.parent/'report_manifest.json',dict(active_package='Co_Directed',manifest='Co_Directed/report_manifest.json',plume_frame='x_p,y_p',obsolete_prediction_displays_excluded=True))
    (out.parent/'README.md').write_text('Active prediction: Co_Directed/README.md, using full-support plume-aligned measured PZT fields. Reproduce from project root: .\\.venv\\Scripts\\python.exe -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook\n',encoding='utf-8')
    return result
