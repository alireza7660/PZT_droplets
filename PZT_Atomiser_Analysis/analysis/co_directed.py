"""Upward co-directed prescribed-mean-field prediction; no source data modified."""
import copy
import json
from pathlib import Path
import numpy as np
from .combined_flow_prediction import ROOT,PROJECT,sources_and_arrays,write_csv,write_json,sha
from .vertical_profile import extract_vertical_profile
from .droplet_dynamics import relaxation_time,acceleration
from .prescribed_gas import MeasuredGasField,integrate_droplet,trajectory_quantities,UnsupportedPosition
from .plotting import plt


def export(fig,out,stem):
    for ext in ('pdf','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)


def context():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());pzt,jet,xj,yj,D0,sources=sources_and_arrays(cfg)
    parent=ROOT/cfg['output'];out=parent/'Co_Directed'
    for folder in ('data','figures','validation'):(out/folder).mkdir(parents=True,exist_ok=True)
    props=json.loads((parent/'fluid_properties.json').read_text())
    profile,peak=extract_vertical_profile(ROOT/cfg['pzt_statistics']);col=peak['native_column']
    paired=np.isfinite(pzt['mean_u'])&np.isfinite(pzt['mean_v'])
    low=float(np.nanmin(pzt['y_mm'][paired]));valid_gas=np.isfinite(jet['U'])&np.isfinite(jet['V'])
    first=int(np.flatnonzero(np.any(valid_gas,axis=0))[0])
    offset=cfg['co_directed']['relative_axial_offset_mm']
    origin=np.array([peak['x_mm']*1e-3,low*1e-3-xj[first]+offset*1e-3])
    field=MeasuredGasField(xj,yj,jet['U'],jet['V'],origin,90)
    alignment=dict(rotation_degrees=90,formula_positions='x=x_reference-Yj; y=Xj+y_PZT_low-Xj_first+axial_offset',
        formula_velocities='u_g=-V_ZNMF; v_g=+U_ZNMF',plume_reference_x_mm=peak['x_mm'],
        PZT_reference_native_column=col,PZT_lowest_plane_y_mm=low,ZNMF_first_supported_Xj_mm=float(xj[first]*1e3),
        relative_axial_offset_mm=offset,translation_m=origin.tolist(),
        interpretation='Measurement-plane alignment assumption, not coincident physical exits; nominal axes aligned, measured centreline wandering retained')
    return cfg,pzt,jet,D0,sources,props,out,profile,peak,field,alignment


def run():
    cfg,pzt,jet,D0,sources,props,out,profile,peak,field,alignment=context();controls=copy.deepcopy(cfg)
    controls['solver']['maximum_time_s']=cfg['co_directed']['maximum_time_s']
    controls['solver']['time_step_s']=cfg['co_directed'].get('time_step_s',cfg['solver']['time_step_s'])
    d=float(cfg['droplet_diameters_um'][0])*1e-6;tau=float(relaxation_time(d,props))
    if len(cfg['droplet_diameters_um'])!=1:raise ValueError('Specify one nominal diameter for this compact co-directed package')
    protected={s['path']:s['sha256'] for s in sources['sources']}
    final=(ROOT/cfg['pzt_statistics']).parent.parent;parent=out.parent
    for f in final.rglob('*'):
        if f.is_file() and parent not in f.resolve().parents and f.suffix.lower() in ('.npz','.csv','.png','.pdf','.svg'):protected[str(f.resolve())]=sha(f)
    columns=[peak['native_column']+int(o) for o in cfg['co_directed']['release_column_offsets']]
    # Lowest common native PZT plane for which every selected release has gas support.
    release_row=None
    for row in np.argsort(pzt['y_mm'][:,peak['native_column']]):
        try:
            for j in columns:
                if not np.isfinite([pzt['mean_u'][row,j],pzt['mean_v'][row,j]]).all():raise ValueError('missing PZT pair')
                field.sample([pzt['x_mm'][row,j]*1e-3,pzt['y_mm'][row,j]*1e-3])
        except (ValueError,UnsupportedPosition):continue
        release_row=int(row);break
    if release_row is None:raise ValueError('No common supported measured PZT release plane for configured axial offset')
    records=[];histories=[];checks=[];force_intervals=[]
    colors=plt.cm.viridis(np.linspace(.05,.85,len(columns)))
    for number,(j,color) in enumerate(zip(columns,colors),1):
        initial=np.array([pzt['x_mm'][release_row,j]*1e-3,pzt['y_mm'][release_row,j]*1e-3,pzt['mean_u'][release_row,j],pzt['mean_v'][release_row,j]])
        tr,reason=integrate_droplet(field,initial,d,props,controls)
        fine,fine_reason=integrate_droplet(field,initial,d,props,controls,.5)
        if len(tr)<2:raise ValueError('No supported integration')
        end=min(tr[-1,0],fine[-1,0]);t=np.linspace(0,end,201)
        errors=[float(np.max(abs(np.interp(t,tr[:,0],tr[:,k])-np.interp(t,fine[:,0],fine[:,k])))) for k in range(1,5)]
        if max(errors[:2])>cfg['solver']['position_convergence_tolerance_m'] or max(errors[2:])>cfg['solver']['velocity_convergence_tolerance_m_s']:raise AssertionError(errors)
        q=trajectory_quantities(field,tr,d,props,D0,cfg)
        a=acceleration(tr[:,3:5],d,np.c_[q['gas_u_m_s'],q['gas_v_m_s']],props,cfg)
        g_eff=props['g_m_s2']*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3'])
        q['net_axial_acceleration_m_s2']=a[:,1];q['drag_axial_acceleration_m_s2']=a[:,1]+g_eff
        q['drag_upward']=q['drag_axial_acceleration_m_s2']>0;q['net_acceleration_upward']=a[:,1]>0
        q['Stk_D']=q['St_D_local'];q['horizontal_displacement_mm']=q['x_mm']-q['x_mm'][0]
        q['vertical_displacement_mm']=q['y_mm']-q['y_mm'][0]
        first_down=np.flatnonzero(q['droplet_v_m_s']<=0)
        end_up=int(first_down[0])+1 if len(first_down) else len(tr)
        peak_index=int(np.argmax(q['y_mm'][:end_up]));ident=f'release_{number:02d}'
        write_csv(out/'data'/(ident+'.csv'),q)
        label=rf'$x_0={initial[0]*1e3:.3f}\,\mathrm{{mm}}$'
        record=dict(release=ident,row=release_row,column=j,x0_mm=float(initial[0]*1e3),y0_mm=float(initial[1]*1e3),
            u0_m_s=float(initial[2]),v0_m_s=float(initial[3]),initial_gas_u_m_s=float(q['gas_u_m_s'][0]),
            initial_gas_v_m_s=float(q['gas_v_m_s'][0]),max_height_mm=float(q['y_mm'][peak_index]),
            time_at_max_height_s=float(q['t_s'][peak_index]),end_time_s=float(tr[-1,0]),termination=reason,
            end_v_m_s=float(q['droplet_v_m_s'][-1]),maximum_particle_Re=float(q['Re_in_plane'].max()),
            initial_net_axial_acceleration_m_s2=float(a[0,1]),any_upward_drag=bool(q['drag_upward'].any()),any_net_upward_acceleration=bool(q['net_acceleration_upward'].any()))
        records.append(record);histories.append(dict(id=ident,q=q,color=color,label=label,peak=peak_index))
        checks.append(dict(release=ident,time_step_s=controls['solver']['time_step_s'],max_errors_x_y_u_v=errors,half_step_termination=fine_reason))
        for definition,mask in [('upward_drag',q['drag_upward']),('net_upward_acceleration',q['net_acceleration_upward']),
                                ('upward_drag_but_net_downward',(q['drag_axial_acceleration_m_s2']>0)&(a[:,1]<0)),
                                ('downward_drag',q['drag_axial_acceleration_m_s2']<0),('net_downward_acceleration',a[:,1]<0)]:
            edges=np.diff(np.r_[False,mask,False].astype(int));starts=np.flatnonzero(edges==1);stops=np.flatnonzero(edges==-1)-1
            for start,stop in zip(starts,stops):force_intervals.append(dict(release=ident,condition=definition,start_time_s=float(q['t_s'][start]),end_time_s=float(q['t_s'][stop]),min_y_mm=float(q['y_mm'][start:stop+1].min()),max_y_mm=float(q['y_mm'][start:stop+1].max())))
    write_csv(out/'data/trajectory_summary.csv',{k:[r[k] for r in records] for k in records[0]})
    if force_intervals:write_csv(out/'data/force_direction_intervals.csv',{k:[r[k] for r in force_intervals] for k in force_intervals[0]})
    lower=max(float(h['q']['y_mm'][0]) for h in histories);upper=min(float(h['q']['y_mm'][h['peak']]) for h in histories)
    if upper<=lower:raise ValueError('No common ascending supported height interval')
    comparisons=[]
    for h in histories:
        q=h['q'];n=h['peak']+1;y=q['y_mm'][:n]
        values={key:[float(np.interp(target,y,q[key][:n])) for target in (lower,upper)] for key in ('t_s','droplet_v_m_s','gas_v_m_s')}
        comparisons.append(dict(release=h['id'],common_lower_y_mm=lower,common_upper_y_mm=upper,
            travel_time_s=values['t_s'][1]-values['t_s'][0],initial_v_m_s=values['droplet_v_m_s'][0],
            final_v_m_s=values['droplet_v_m_s'][1],delta_v_m_s=values['droplet_v_m_s'][1]-values['droplet_v_m_s'][0],
            gas_v_at_lower_m_s=values['gas_v_m_s'][0],gas_v_at_upper_m_s=values['gas_v_m_s'][1]))
    write_csv(out/'data/common_height_comparison.csv',{k:[r[k] for r in comparisons] for k in comparisons[0]})
    write_csv(out/'data/PZT_only_fixed_line_profile.csv',profile)
    stems=plots(out,field,histories,profile,peak,g_eff)
    from .co_directed_report import dimensionless,write_report
    numbers,jet_stem=dimensionless(out,cfg,props,D0,alignment);stems.append(jet_stem)
    result=dict(successfully_integrated_cases=len(records),alignment=alignment,release_row=release_row,release_columns=columns,
        common_interval_y_mm=[lower,upper],common_interval_extent_mm=upper-lower,
        common_travel_time_range_ms=[min(c['travel_time_s'] for c in comparisons)*1e3,max(c['travel_time_s'] for c in comparisons)*1e3],
        common_delta_v_range_m_s=[min(c['delta_v_m_s'] for c in comparisons),max(c['delta_v_m_s'] for c in comparisons)],
        peak_height_range_mm=[min(r['max_height_mm'] for r in records),max(r['max_height_mm'] for r in records)],
        termination_reasons={reason:sum(r['termination']==reason for r in records) for reason in {r['termination'] for r in records}},
        initial_gas_axial_velocity_range_m_s=[min(r['initial_gas_v_m_s'] for r in records),max(r['initial_gas_v_m_s'] for r in records)],
        reversed_cases=sum(r['end_v_m_s']<0 for r in records),dimensionless=numbers)
    result['cases_with_upward_drag']=sum(r['any_upward_drag'] for r in records)
    result['cases_with_net_upward_acceleration']=sum(r['any_net_upward_acceleration'] for r in records)
    write_json(out/'results_summary.json',result);write_json(out/'configuration_used.json',cfg)
    write_json(out/'source_manifest.json',dict(**sources,co_directed_alignment=alignment,saved_properties=props,
        assumptions='Nominal axes co-directed upward; signed local gas velocities retained, not forced positive'))
    write_json(out/'validation/convergence.json',checks);write_json(out/'validation/authoritative_hashes.json',protected)
    write_json(out/'report_manifest.json',dict(figures=stems,scalar_table='data/scalar_quantities.csv',
        numerical_data=sorted(p.name for p in (out/'data').glob('*.csv')),text=['methods_results.tex','captions.tex'],archive_manifest=cfg['co_directed']['archive_manifest']))
    write_report(out,cfg,result,stems)
    for path,digest in protected.items():assert sha(path)==digest,path
    write_json(parent/'prediction_status.json',dict(status='co_directed_prediction_complete',prediction_case_count=len(records),active_package='Co_Directed',combined_experiment_conducted=False))
    write_json(parent/'report_manifest.json',dict(active_package='Co_Directed',manifest='Co_Directed/report_manifest.json',archived_side_jet_excluded=True))
    (parent/'README.md').write_text('Active package: Co_Directed/. Horizontal side-jet outputs are archived. Regenerate using .\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook from the workspace root.\n',encoding='utf-8')
    return result


def _branches(ax,h,key,gas=False):
    q=h['q'];k=h['peak'];color=h['color'];style=':' if gas else '-'
    ax.plot(q['y_mm'][:k+1],q[key][:k+1],color=color,ls=style,lw=1.2,label=h['label'] if not gas else None)
    if k+1<len(q['y_mm']):ax.plot(q['y_mm'][k:],q[key][k:],color=color,ls='-.' if gas else '--',lw=1.1)


def plots(out,field,histories,profile,peak,g_eff):
    from matplotlib.colors import Normalize
    xx,yy,ug,vg=field.common_mesh();extent=float(np.nanmax(abs(vg)));stems=[]
    x_all=np.concatenate([h['q']['x_mm'] for h in histories]);y_all=np.concatenate([h['q']['y_mm'] for h in histories])
    for stem,detail in [('trajectories_upward_jet',False),('trajectories_entry_detail',True)]:
        fig,ax=plt.subplots(figsize=(6.4,5.2) if not detail else (8,3.6),layout='constrained')
        m=ax.contourf(xx*1e3,yy*1e3,np.ma.masked_invalid(vg),levels=np.linspace(-extent,extent,41),cmap='RdBu_r',norm=Normalize(-extent,extent),corner_mask=False)
        if detail:fig.colorbar(m,ax=ax,label=r'$\overline{v}_g\ (\mathrm{m\,s^{-1}})$',location='top',orientation='horizontal',aspect=50,shrink=.8,pad=.05)
        else:fig.colorbar(m,ax=ax,label=r'$\overline{v}_g\ (\mathrm{m\,s^{-1}})$',shrink=.8)
        for h in histories:
            q=h['q'];ax.plot(q['x_mm'],q['y_mm'],color=h['color'],label=h['label'],lw=1.4)
            ax.plot(q['x_mm'][0],q['y_mm'][0],'o',ms=3,color=h['color'])
        if detail:ax.set(xlim=(x_all.min()-.12,x_all.max()+.12),ylim=(y_all.min()-.015,y_all.max()+.025))
        ax.set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$');ax.set_aspect('equal')
        if detail:fig.legend(*ax.get_legend_handles_labels(),frameon=False,fontsize=8,ncol=4,loc='lower center',bbox_to_anchor=(.5,-.025))
        export(fig,out,stem);stems.append(stem)
    fig,ax=plt.subplots(figsize=(7.2,4.6),layout='constrained')
    for h in histories:_branches(ax,h,'droplet_v_m_s');_branches(ax,h,'gas_v_m_s',True)
    ax.axhline(0,color='0.6',lw=.6);ax.set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$v_d,\ v_g\ (\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=7,ncol=2)
    export(fig,out,'axial_velocity_along_trajectories');stems.append('axial_velocity_along_trajectories')
    fig,axes=plt.subplots(2,1,figsize=(7.2,6),sharex=True,layout='constrained')
    for h in histories:
        _branches(axes[0],h,'slip_v_m_s');_branches(axes[1],h,'net_axial_acceleration_m_s2')
        # Drag alone is plotted as faint dotted traces to distinguish it from net acceleration.
        q=h['q'];axes[1].plot(q['y_mm'],q['drag_axial_acceleration_m_s2'],':',color=h['color'],alpha=.6,lw=.8)
    for ax in axes:ax.axhline(0,color='0.5',lw=.7)
    axes[1].axhline(g_eff,color='black',ls=':',lw=.7)
    axes[0].set_ylabel(r'$v_g-v_d\ (\mathrm{m\,s^{-1}})$');axes[0].legend(frameon=False,fontsize=7,ncol=2)
    axes[1].set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$a_y,\ a_{D,y}\ (\mathrm{m\,s^{-2}})$',yscale='symlog');axes[1].set_yscale('symlog',linthresh=20)
    export(fig,out,'axial_slip_and_acceleration');stems.append('axial_slip_and_acceleration')
    fig,axes=plt.subplots(2,1,figsize=(7.2,5.6),sharex=True,layout='constrained')
    for h in histories:_branches(axes[0],h,'Re_in_plane');_branches(axes[1],h,'Stk_D')
    axes[0].set_ylabel(r'$Re_p\ (-)$');axes[0].legend(frameon=False,fontsize=7,ncol=2)
    axes[1].set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$Stk_D\ (-)$')
    export(fig,out,'particle_Re_and_local_Stokes');stems.append('particle_Re_and_local_Stokes')
    fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
    ax.plot(profile['y_mm'],profile['mean_v'],color='#d62728',label=rf'Measured fixed $x={peak["x_mm"]:.5f}$ mm')
    ax.set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$\overline{v}_d\ (\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=9)
    export(fig,out,'PZT_only_measured_axial_profile');stems.append('PZT_only_measured_axial_profile')
    return stems


def notebook_cells():
    import nbformat as nb
    first=nb.v4.new_markdown_cell('## Upward co-directed prediction\n\nSeparate PZT droplet and ZNMF gas measurements; nominal axes aligned using explicitly assumed coincident measurement planes, not physical exits. Signed local gas reversals and measured centreline wandering are retained. Side-jet displays are removed.')
    code=nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.co_directed import run as run_codirected
prediction=run_codirected()
print(json.dumps(prediction,indent=2))
package=output/'Combined_Flow_Prediction/Co_Directed'
for name in ['trajectories_upward_jet','trajectories_entry_detail','axial_velocity_along_trajectories']:
    display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    text=nb.v4.new_markdown_cell('## Forces, common-interval comparison and supporting dimensionless numbers\n\nSolid/dashed particle curves distinguish initial ascent and subsequent descent; gas traces are dotted/dash-dot. The red PZT-only ensemble baseline is shown separately. Scalar response/frequency values are tabulated, not plotted as constant curves.')
    figures=nb.v4.new_code_cell('''print((package/'data/common_height_comparison.csv').read_text())
print((package/'data/scalar_quantities.csv').read_text())
for name in ['axial_slip_and_acceleration','particle_Re_and_local_Stokes','PZT_only_measured_axial_profile','jet_local_dimensionless_numbers']:
    display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    for c in (first,code,text,figures):c.metadata['co_directed']=True
    first.metadata['combined_flow_prediction']=True;code.metadata['combined_flow_prediction']=True
    return [first,code,text,figures]
