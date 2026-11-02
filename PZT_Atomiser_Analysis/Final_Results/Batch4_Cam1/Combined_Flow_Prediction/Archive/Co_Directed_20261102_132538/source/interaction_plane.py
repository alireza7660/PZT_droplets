"""Trace the archived comparison to authoritative profiles and build release audit."""
import json
import numpy as np
from .combined_flow_prediction import ROOT,PROJECT,write_csv,write_json,sha
from .prescribed_gas import MeasuredGasField
from .droplet_dynamics import acceleration,reynolds
from .plotting import plt

COMPARISON=PROJECT/'Archive/20261102T095744_AEDT_superseded_interaction/Jet_PZT_Comparison'


def bracket(x,point,valid):
    i=int(np.searchsorted(x,point,side='right'))-1;i=min(i,len(x)-2)
    if i<0 or not (x[i]<=point<=x[i+1]) or not valid[i:i+2].all():raise ValueError('Unsupported interpolation interval')
    w=float((point-x[i])/(x[i+1]-x[i]));return i,w


def supported_value(x,values,valid,point):
    i,w=bracket(x,point,valid);return float((1-w)*values[i]+w*values[i+1])


def prepare(cfg,pzt,jet,xj,yj,D0,profile,peak,props,out,write=True):
    gas_cfg=(ROOT/cfg['znmf_configuration']).resolve()
    gas_profile=gas_cfg.parent.parent/'tables/FOV_1_2_3_centreline_and_width.csv'
    pzt_profile=(ROOT/cfg['pzt_statistics']).parent/'mean_v_vertical_profile.csv'
    z=np.genfromtxt(gas_profile,names=True,delimiter=',',dtype=None,encoding='utf-8')
    original_p=np.genfromtxt(pzt_profile,names=True,delimiter=',')
    for key in ('y_mm','mean_v','x_mm','native_row','native_column'):np.testing.assert_array_equal(original_p[key],profile[key])
    diameter_mm=float(json.loads(gas_cfg.read_text())['diameter_mm'])
    sx=z['x_over_D0']*diameter_mm;uc=z['Uc_m_per_s'];py=profile['y_mm'];pv=profile['mean_v']
    valid=(z['status']=='measured')&np.isfinite(uc)&(uc>0)&np.isfinite(z['yc_absolute_over_D0'])
    pvalid=np.isfinite(pv);offset=float(cfg['co_directed']['relative_axial_offset_mm'])
    # Original figure used each independent measured coordinate on 0..27 mm.
    zg=(z['status']=='measured')&np.isfinite(uc)&(sx>=0)&(sx<=27)
    pg=np.isfinite(pv)&(py>=0)&(py<=27)
    archived_z=np.genfromtxt(COMPARISON/'znmf_plotted_data.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    archived_p=np.genfromtxt(COMPARISON/'pzt_plotted_data.csv',names=True,delimiter=',')
    np.testing.assert_array_equal(sx[zg],archived_z['s_mm']);np.testing.assert_array_equal(uc[zg],archived_z['mean_axial_velocity_m_s'])
    np.testing.assert_array_equal(py[pg],archived_p['s_mm']);np.testing.assert_array_equal(pv[pg],archived_p['mean_axial_velocity_m_s'])
    nodes=np.unique(np.r_[sx[(sx>=0)&(sx<=27)],py[(py-offset>=0)&(py-offset<=27)]-offset])
    roots=[]
    for left,right in zip(nodes[:-1],nodes[1:]):
        mid=(left+right)/2
        try:
            gi,_=bracket(sx,mid,valid);pi,_=bracket(py,mid+offset,pvalid)
            g=np.interp([left,right],sx[gi:gi+2],uc[gi:gi+2]);p=np.interp(np.array([left,right])+offset,py[pi:pi+2],pv[pi:pi+2])
        except ValueError:continue
        diff=g-p
        if diff[0]<=0<=diff[1] and diff[1]>diff[0]:
            root=float(left-diff[0]*(right-left)/(diff[1]-diff[0]))
            if roots and abs(roots[-1]['ZNMF_downstream_mm']-root)<1e-10:continue
            roots.append(dict(ZNMF_downstream_mm=root,PZT_height_mm=root+offset,
                common_velocity_m_s=supported_value(sx,uc,valid,root),gas_lower_source_column=gi,
                PZT_lower_profile_index=pi,bracket_lower_mm=float(left),bracket_upper_mm=float(right)))
    if not roots:raise ValueError('No supported gas-overtaking crossing in the actual comparison domain')
    crossing=roots[0];sg=crossing['ZNMF_downstream_mm'];height=crossing['PZT_height_mm']
    yc=supported_value(sx,z['yc_absolute_over_D0'],valid,sg)
    source_y0=json.loads((ROOT/cfg['znmf_nozzle_reference']).read_text())['geometric_nozzle_centre_registered_D0']
    local_yc=(yc-source_y0)*D0
    origin=np.array([peak['x_mm']*1e-3+local_yc,(height-sg)*1e-3])
    field=MeasuredGasField(xj,yj,jet['U'],jet['V'],origin,90)
    alignment=dict(rotation_degrees=90,formula_positions='x=x_PZT_ref + Yj_c_at_interaction - Yj; y=y_PZT_interaction + Xj - s_ZNMF_interaction',
        formula_velocities='u_g=-V_ZNMF; v_g=+U_ZNMF',plume_reference_x_mm=peak['x_mm'],PZT_reference_native_column=peak['native_column'],
        PZT_interaction_height_mm=height,ZNMF_interaction_downstream_mm=sg,ZNMF_measured_centreline_y_registered_over_D0=yc,
        ZNMF_interaction_centreline_Yj_mm=local_yc*1e3,relative_axial_offset_mm=offset,translation_m=origin.tolist(),
        interpretation='Independent source profile coordinates paired as in archived comparison; measured centrelines aligned only at interaction plane, not physical exits; wandering and transverse offsets retained')
    states=[];initial_table=[];stencils=[];d=cfg['droplet_diameters_um'][0]*1e-6
    for j in [peak['native_column']+int(o) for o in cfg['co_directed']['release_column_offsets']]:
        rows=np.argsort(pzt['y_mm'][:,j]);y=pzt['y_mm'][rows,j]
        vvalid=np.isfinite(pzt['mean_u'][rows,j])&np.isfinite(pzt['mean_v'][rows,j])
        k,w=bracket(y,height,vvalid);r0,r1=int(rows[k]),int(rows[k+1])
        u=float((1-w)*pzt['mean_u'][r0,j]+w*pzt['mean_u'][r1,j]);v=float((1-w)*pzt['mean_v'][r0,j]+w*pzt['mean_v'][r1,j])
        initial=np.array([pzt['x_mm'][r0,j]*1e-3,height*1e-3,u,v]);gas=field.sample(initial[:2])
        local=field.to_local(initial[:2]);ix,wx=bracket(xj,local[0],np.ones(len(xj),bool));iy,wy=bracket(yj,local[1],np.ones(len(yj),bool))
        acc=acceleration(initial[2:],d,gas,props,cfg);g_eff=props['g_m_s2']*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3'])
        row=dict(PZT_native_column=j,PZT_lower_native_row=r0,PZT_upper_native_row=r1,PZT_upper_weight=w,
            PZT_lower_y_mm=float(y[k]),PZT_upper_y_mm=float(y[k+1]),PZT_source_x_mm=float(initial[0]*1e3),
            PZT_source_height_mm=height,ZNMF_source_downstream_mm=float(local[0]*1e3),ZNMF_source_Yj_mm=float(local[1]*1e3),
            ZNMF_source_x_over_D0=float(local[0]/D0),ZNMF_source_registered_y_over_D0=float(local[1]/D0+source_y0),
            gas_lower_grid_column=ix,gas_lower_grid_row=iy,gas_streamwise_upper_weight=wx,gas_transverse_upper_weight=wy,
            common_x_mm=float(initial[0]*1e3),common_y_mm=height,measured_droplet_u_m_s=u,measured_droplet_v_m_s=v,
            sampled_gas_u_m_s=float(gas[0]),sampled_gas_v_m_s=float(gas[1]),plotted_centreline_Uc_m_s=crossing['common_velocity_m_s'],
            sampled_minus_plotted_axial_m_s=float(gas[1]-crossing['common_velocity_m_s']),slip_u_m_s=float(gas[0]-u),slip_v_m_s=float(gas[1]-v),
            initial_Re_p=float(reynolds(d,initial[2:],gas,props)),initial_drag_axial_acceleration_m_s2=float(acc[1]+g_eff),
            initial_net_axial_acceleration_m_s2=float(acc[1]))
        initial_table.append(row);states.append(dict(column=j,initial=initial,PZT_lower_row=r0,PZT_upper_row=r1,upper_weight=w))
        stencils.append(dict(column=j,PZT_lower_row=r0,PZT_upper_row=r1,PZT_upper_weight=w,
            PZT_u_source=[float(pzt['mean_u'][r0,j]),float(pzt['mean_u'][r1,j])],PZT_v_source=[float(pzt['mean_v'][r0,j]),float(pzt['mean_v'][r1,j])],
            gas_U_corners=jet['U'][iy:iy+2,ix:ix+2].tolist(),gas_V_corners=jet['V'][iy:iy+2,ix:ix+2].tolist(),
            gas_corner_weight_matrix=[[(1-wx)*(1-wy),wx*(1-wy)],[(1-wx)*wy,wx*wy]]))
    centre=next(r for r in initial_table if r['PZT_native_column']==peak['native_column'])
    np.testing.assert_allclose(centre['measured_droplet_v_m_s'],crossing['common_velocity_m_s'],atol=1e-12)
    evidence=dict(selected_crossing=crossing,all_supported_gas_overtaking_crossings=roots,
        comparison_reproduction=dict(gas_points=int(zg.sum()),PZT_points=int(pg.sum()),source_curve_values_exactly_match_archived_CSVs=True,
            archived_script=str(COMPARISON/'plot_comparison.py'),archived_readme=str(COMPARISON/'README.md')),
        source_paths=dict(PZT_profile=str(pzt_profile),ZNMF_centreline=str(gas_profile)),
        processing_evidence=dict(ZNMF_configuration=str(gas_cfg),ZNMF_centreline_extraction_code=str(gas_cfg.parent/'pipeline.py'),
            ZNMF_extraction_code_sha256=sha(gas_cfg.parent/'pipeline.py'),PZT_profile_code=str(PROJECT/'analysis/vertical_profile.py')),
        centre_release_sampled_minus_plotted_Uc_m_s=centre['sampled_minus_plotted_axial_m_s'],
        reason_for_gas_difference='Uc is evaluated at each native column measured centreline, then interpolated downstream. Field samples bilinearly interpolate U,V at one interpolated centreline/transverse position. These operations differ when the centreline wanders. Off-axis releases sample different transverse profiles. See ZNMF pipeline.extract_centreline: quadratic position estimate followed by linear velocity interpolation.',
        initial_conditions_exported_before_integration=True,calibration_not_reapplied=True)
    if write:
        write_csv(out/'data/initial_conditions.csv',{k:[r[k] for r in initial_table] for k in initial_table[0]})
        write_json(out/'data/initial_interpolation_stencils.json',stencils);write_json(out/'interaction_plane_audit.json',evidence)
        write_csv(out/'data/comparison_source_ZNMF.csv',dict(ZNMF_downstream_mm=sx[zg],Uc_m_s=uc[zg]))
        write_csv(out/'data/comparison_source_PZT.csv',dict(PZT_height_mm=py[pg],mean_v_m_s=pv[pg]))
        fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
        ax.plot(sx[zg],uc[zg],color='black',ls='--',label='Measured ZNMF centreline')
        ax.plot(py[pg],pv[pg],color='#d62728',label='Measured PZT-only fixed line')
        ax.plot(sg,crossing['common_velocity_m_s'],'o',color='black',ms=4,label=rf'Crossing $s={sg:.5f}$ mm')
        ax.set(xlim=(0,27),ylim=(0,2.4),xlabel=r'$s\ (\mathrm{mm})$',ylabel=r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
        for ext in ('pdf','png'):fig.savefig(out/'figures'/('source_comparison_reproduced.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig)
    return field,alignment,states,evidence
