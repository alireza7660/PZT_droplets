"""Rigid jet placement beneath the first supported measured plume-axis plane."""
import json
import numpy as np
from .combined_flow_prediction import ROOT,write_json,write_csv
from .prescribed_gas import MeasuredGasField,UnsupportedPosition
from .interaction_plane import bracket
from .droplet_dynamics import acceleration,reynolds,relaxation_time,drag_factor


def prepare_base(ctx,write=True):
    from .aligned_codirected import pzt_sample
    cfg=ctx['cfg'];p=ctx['profile'];base=ctx['base'];target=float(cfg['co_directed']['entrance_centreline_velocity_m_s']);sx=ctx['sx'];uc=ctx['uc'];valid=ctx['valid']
    peak=int(np.nanargmax(np.where(valid,uc,np.nan)));brackets=[]
    for i in range(peak):
        if valid[i:i+2].all() and uc[i]<=target<=uc[i+1] and uc[i+1]>uc[i]:
            w=float((target-uc[i])/(uc[i+1]-uc[i]));s=float(sx[i]+w*(sx[i+1]-sx[i]));brackets.append(dict(lower_native_column=i,upper_native_column=i+1,upper_weight=w,x_j0_mm=s,lower_Xj_mm=float(sx[i]),upper_Xj_mm=float(sx[i+1]),lower_Uc_m_s=float(uc[i]),upper_Uc_m_s=float(uc[i+1])))
    if not brackets:raise ValueError('No adjacent valid developing-branch centreline samples bracket target')
    selected=brackets[0];i=selected['lower_native_column'];w=selected['upper_weight'];s0=selected['x_j0_mm']
    curve_value=float((1-w)*uc[i]+w*uc[i+1]);np.testing.assert_allclose(curve_value,target,rtol=0,atol=1e-12)
    nozzle_y=ctx['sources']['znmf']['nozzle_y_registered_D0'];yc_registered=float((1-w)*ctx['gas']['yc_absolute_over_D0'][i]+w*ctx['gas']['yc_absolute_over_D0'][i+1]);yc=(yc_registered-nozzle_y)*ctx['D0']
    field=MeasuredGasField(ctx['xj'],ctx['yj'],ctx['jet']['U'],ctx['jet']['V'],[yc,-s0*1e-3],90.)
    drop,reference_stencil=pzt_sample(ctx,[0,base['y_p_mm']]);gas=field.sample([0,0])
    entry=dict(ZNMF_downstream_mm=s0,ZNMF_local_centreline_Y_mm=float(yc*1e3),measured_Uc_m_s=curve_value,target_centreline_Uc_m_s=target,sampled_gas_v_p_m_s=float(gas[1]),sampled_gas_u_p_m_s=float(gas[0]),sampled_minus_curve_v_m_s=float(gas[1]-curve_value),PZT_lab_x_mm=base['lab_x_mm'],PZT_lab_y_mm=base['lab_y_mm'],x_p_mm=0.,y_p_mm=base['y_p_mm'],h_mm=0.,PZT_axis_u_p_m_s=float(drop[0]),PZT_axis_v_p_m_s=float(drop[1]),PZT_saved_axis_profile_interpolated_v_p_m_s=float(p['mean_v_p_m_s'][base['profile_index']]),PZT_axis_native_field_minus_profile_m_s=float(drop[1]-p['mean_v_p_m_s'][base['profile_index']]))
    np.testing.assert_allclose(field.to_local([0,0]),[s0*1e-3,yc],atol=1e-14)
    gravity=np.array(cfg['gravity_vector_m_s2'])*(1-ctx['props']['rho_air_kg_m3']/ctx['props']['rho_water_kg_m3']);d=cfg['droplet_diameters_um'][0]*1e-6
    spacing=float(np.median(np.diff(ctx['a']['x_p_mm'],axis=1)));states=[];table=[];stencils=[];candidates=[]
    for off in cfg['co_directed']['release_column_offsets']:
        xp=float(off*spacing);point_plume=np.array([xp,base['y_p_mm']]);lab=point_plume@ctx['Q']+ctx['origin'];candidate=dict(transverse_grid_offset=int(off),x_p_mm=xp,h_mm=0.,y_p_mm=base['y_p_mm'],PZT_lab_x_mm=float(lab[0]),PZT_lab_y_mm=float(lab[1]))
        try:v,stencil=pzt_sample(ctx,point_plume)
        except ValueError as problem:
            candidates.append(dict(**candidate,accepted=False,reason='unsupported_PZT_native_cell: '+str(problem)));continue
        try:g=field.sample([xp*1e-3,0])
        except UnsupportedPosition as problem:
            candidates.append(dict(**candidate,accepted=False,reason='unsupported_gas_native_cell: '+str(problem)));continue
        number=len(states)+1;initial=np.r_[xp*1e-3,0.,v];net=acceleration(v,d,g,ctx['props'],cfg);drag=net-gravity
        local=field.to_local(initial[:2]);j,wx=bracket(ctx['xj'],local[0],np.ones(len(ctx['xj']),bool));k,wy=bracket(ctx['yj'],local[1],np.ones(len(ctx['yj']),bool));weights=np.array([[(1-wx)*(1-wy),wx*(1-wy)],[(1-wx)*wy,wx*wy]])
        record=dict(release_number=number,is_reference=off==0,**candidate,ZNMF_downstream_mm=float(local[0]*1e3),ZNMF_local_Y_mm=float(local[1]*1e3),PZT_native_row_0=stencil['native_rows_in_weight_order'][0],PZT_native_row_1=stencil['native_rows_in_weight_order'][1],PZT_native_col_0=stencil['native_columns'][0],PZT_native_col_1=stencil['native_columns'][1],PZT_weight_00=float(stencil['weights'][0][0]),PZT_weight_01=float(stencil['weights'][0][1]),PZT_weight_10=float(stencil['weights'][1][0]),PZT_weight_11=float(stencil['weights'][1][1]),droplet_u_p_m_s=float(v[0]),droplet_v_p_m_s=float(v[1]),gas_u_p_m_s=float(g[0]),gas_v_p_m_s=float(g[1]),measured_centreline_Uc_m_s=curve_value,sampled_minus_centreline_v_m_s=float(g[1]-curve_value),slip_u_p_m_s=float(g[0]-v[0]),slip_v_p_m_s=float(g[1]-v[1]),initial_a_x_p_m_s2=float(net[0]),initial_a_y_p_m_s2=float(net[1]),initial_drag_x_p_m_s2=float(drag[0]),initial_drag_y_p_m_s2=float(drag[1]),Re_p=float(reynolds(d,v,g,ctx['props'])))
        stencil.update(release_number=number,point_common_x_p_h_mm=[xp,0.],gas_lower_native_row=k,gas_lower_native_column=j,gas_weights=weights.tolist(),gas_U_corners=ctx['jet']['U'][k:k+2,j:j+2].tolist(),gas_V_corners=ctx['jet']['V'][k:k+2,j:j+2].tolist())
        states.append(dict(initial=initial,reference=off==0,number=number));table.append(record);stencils.append(stencil);candidates.append(dict(**candidate,accepted=True,reason='finite_PZT_and_gas_native_cells'))
    if not any(s['reference'] for s in states):raise ValueError('Axis reference cannot be initialised at first supported measurement plane')
    alignment=dict(entry=entry,base_measurement=base,developing_branch_bracket=selected,gas_to_plume_rotation_matrix=field.R.tolist(),gas_nozzle_translation_plume_m=field.origin.tolist(),PZT_lab_to_plume_Q=ctx['Q'].tolist(),PZT_reference_lab_mm=ctx['origin'].tolist(),gravity_plume_m_s2=cfg['gravity_vector_m_s2'],rigid_placement='[x_p,h]=R90*[Xj,Yj]+[Yc(x_j0),-x_j0]; [u_g,p,v_g,p]=[-Vj,Uj]; h=y_p-y_p_base. No stretching or downstream recentering.',gas_streamwise_to_plume_formula='Xj_mm=x_j0_mm+h_mm',inverse_to_original_plume='[x_p,y_p]=[x_p,h+y_p_base]',inverse_to_lab='[x_lab,y_lab]=Q^T*[x_p,h+y_p_base]+PZT_reference_lab',interpretation='Jet-under-plume prescribed design placement at the first measured axis plane; not the physical atomiser surface or an equating of experimental origins. Previous profile-crossing placement is superseded.')
    audit=dict(placement_mode='plume_base_target_centreline',selected_entry=entry,base_measurement=base,developing_branch_bracket=selected,all_supported_developing_target_brackets=brackets,release_candidates=candidates,accepted_release_count=len(states),rejected_release_count=len(candidates)-len(states),crossing_used_to_select_entry=False,initial_conditions_exported_before_integration=True,calibration_reapplied=False,report_crop_used=False,gas_centreline_difference='The target 0.96 m/s is exact on the measured centreline curve interpolated between adjacent supported stations. Native bilinear U,V sampled at the interpolated centreline can differ because downstream interpolation and centreline extraction do not commute. Both values and raw stencils are reported; the gas field is not rescaled or forced uniform.')
    if write:
        out=ctx['out'];write_csv(out/'data/initial_conditions.csv',{key:[r[key] for r in table] for key in table[0]});write_csv(out/'data/release_candidate_audit.csv',{key:[r[key] for r in candidates] for key in candidates[0]});write_json(out/'data/initial_interpolation_stencils.json',stencils);write_json(out/'interaction_plane_audit.json',audit);write_json(out/'transformation.json',alignment)
        write_csv(out/'data/entry_mapping.csv',dict(frame=['ZNMF local nozzle reference','Original PZT laboratory','Saved plume-aligned','Prediction entrance frame'],coordinate_1_mm=[s0,base['lab_x_mm'],0.,0.],coordinate_2_mm=[float(yc*1e3),base['lab_y_mm'],base['y_p_mm'],0.],axis_names=['Xj; Yj','x; y','x_p; y_p','x_p; h']))
        write_csv(out/'data/PZT_corrected_axis_profile.csv',p);write_csv(out/'data/ZNMF_profile_coordinate_mapping.csv',dict(ZNMF_downstream_mm=sx,h_mm=sx-s0,original_plume_y_p_mm=sx-s0+base['y_p_mm'],measured_Uc_m_s=uc,centreline_valid=valid,status=ctx['gas']['status']))
    return field,states,table,alignment,audit
