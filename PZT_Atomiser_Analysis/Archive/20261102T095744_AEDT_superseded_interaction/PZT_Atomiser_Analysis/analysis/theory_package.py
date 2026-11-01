"""Report package: distinguish Eulerian mean profiles from Lagrangian references.

The dynamical equations remain in theory_comparison.py. Every figure is drawn
from the exported CSV values; a manifest records the exact curve provenance.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from . import theory_comparison as model
from .plotting import plt

FIGURES=(
    'velocity_reference_full','velocity_reference_near_origin',
    'reynolds_measured_velocity_assumed_air','reynolds_reference_trajectories',
    'particle_relaxation_time','diameter_sensitivity',
    'Jet_Interaction_Sensitivity/znmf_local_nozzle_stokes')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    return np.genfromtxt(path,names=True,delimiter=',',encoding='utf-8')


def curve(ax,out,csv_name,xkey,ykey,label,*,color,linestyle='-',linewidth=1.2,
          xfactor=1.,yfactor=1.,filter_column=None,filter_value=None):
    data=read_csv(out/csv_name)
    if filter_column is not None:data=data[data[filter_column]==filter_value]
    x=data[xkey]*xfactor;y=data[ykey]*yfactor
    line,=ax.plot(x,y,label=label,color=color,ls=linestyle,lw=linewidth)
    np.testing.assert_array_equal(line.get_xdata(),x)
    np.testing.assert_array_equal(line.get_ydata(),y)
    return dict(csv=csv_name,x_column=xkey,y_column=ykey,x_factor=xfactor,y_factor=yfactor,
                filter_column=filter_column,filter_value=filter_value,label=label,color=color,
                linestyle=linestyle,linewidth=linewidth,point_count=len(x),
                x_sha256=hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest(),
                y_sha256=hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest())


def export(fig,ax,out,stem,curves,manifest,*,velocity=False,category='reference'):
    if velocity:ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2)
    else:ax.legend(frameon=False,fontsize=8)
    assert not ax.get_title()
    target=out/stem;target.parent.mkdir(parents=True,exist_ok=True)
    entry=dict(stem=stem,category=category,curves=curves,xlim=list(ax.get_xlim()),
               ylim=list(ax.get_ylim()),xscale=ax.get_xscale(),yscale=ax.get_yscale(),
               exports=[stem+'.'+ext for ext in ('pdf','svg','png')],png_dpi=600)
    for ext in ('pdf','svg','png'):fig.savefig(target.with_suffix('.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig);manifest.append(entry)


def axis(xlabel,ylabel):
    fig,ax=plt.subplots(figsize=(6.4,4.5),layout='constrained')
    ax.set(xlabel=xlabel,ylabel=ylabel)
    return fig,ax


def endpoint_record(name,diameter,initial,plotted,status,turning,turn_status,y0):
    if turn_status!='turning_point':raise ValueError('Actual stopping point was not reached; increase time limit')
    return dict(reference=name,assumed_diameter_um=diameter,
                actual_turning_height_mm=float(turning[-1,2]*1e3-y0),
                actual_turning_y_mm=float(turning[-1,2]*1e3),
                actual_turning_time_s=float(turning[-1,0]),
                plotted_endpoint_height_mm=float(plotted[-1,2]*1e3-y0),
                plotted_endpoint_y_mm=float(plotted[-1,2]*1e3),
                plotted_endpoint_time_s=float(plotted[-1,0]),
                plotted_endpoint_reason=status,initial_u_m_s=float(initial[2]),initial_v_m_s=float(initial[3]))


def run(config_path=None):
    cfg_path=Path(config_path or model.PROJECT/'config/theory_comparison.json')
    cfg=json.loads(cfg_path.read_text());p=model.fluid_properties(cfg)
    if np.any(np.asarray(cfg['gas_velocity_m_s'])!=0):raise ValueError('Stationary-air references require zero gas velocity')
    if cfg.get('measured_gas_profile_csv'):raise ValueError('This report explicitly separates stationary-air references from measurements')
    out=model.FINAL/'Theory_Comparison';out.mkdir(exist_ok=True)
    jet=out/'Jet_Interaction_Sensitivity';jet.mkdir(exist_ok=True)
    validation_dir=out/'validation';validation_dir.mkdir(exist_ok=True)
    protected={str(f.relative_to(model.FINAL)):digest(f) for f in model.FINAL.rglob('*')
               if f.is_file() and out not in f.parents and f.suffix.lower() in ('.csv','.npz','.pdf','.svg','.png')}
    active=json.loads((model.PROJECT/'config/analysis.json').read_text())
    np.testing.assert_allclose([active['coordinate_scale'],active['velocity_scale'],active['confirmed_pulse_separation_seconds']],
                               [.00535,.214,25e-6],rtol=0,atol=1e-15)
    profile,maximum=model.extract_vertical_profile(model.FINAL/'data/statistics.npz')
    with np.load(model.FINAL/'data/statistics.npz') as data:
        u=data['mean_u'][profile['native_row'],profile['native_column']]
    y=profile['y_mm'];v=profile['mean_v'];w=np.c_[u,v]
    good=np.isfinite(w).all(axis=1)&(profile['N']>=cfg['minimum_initial_sample_count'])
    first=np.flatnonzero(good)[0];y0=float(y[first]);height=y-y0
    initial=np.array([profile['x_mm'][first]*1e-3,y0*1e-3,u[first],v[first]])
    if initial[3]<=0:raise ValueError('Selected initial measured axial velocity is not upward; do not tune it')
    top=float(y.max())*1e-3;hmax=float(height.max())
    model.write_csv(out/'measured_fixed_x_profile.csv',dict(**profile,mean_u=u,height_above_first_measurement_mm=height))
    g=p['g_m_s2'];t_stop=initial[3]/g;h_stop=initial[3]**2/(2*g)*1e3
    h_end=min(hmax,h_stop);ballistic_h=np.linspace(0,h_end,501)
    ballistic_v=np.sqrt(np.maximum(0,initial[3]**2-2*g*ballistic_h*1e-3))
    ballistic_t=(initial[3]-ballistic_v)/g
    model.write_csv(out/'ballistic_reference.csv',dict(y_mm=ballistic_h+y0,
                    height_above_first_measurement_mm=ballistic_h,t_s=ballistic_t,v_m_s=ballistic_v))
    endpoints=[dict(reference='gravity_only',assumed_diameter_um=None,
        actual_turning_height_mm=float(h_stop),actual_turning_y_mm=float(h_stop+y0),
        actual_turning_time_s=float(t_stop),plotted_endpoint_height_mm=float(h_end),
        plotted_endpoint_y_mm=float(h_end+y0),plotted_endpoint_time_s=float(ballistic_t[-1]),
        plotted_endpoint_reason='domain_top' if hmax<h_stop else 'turning_point',
        initial_u_m_s=float(initial[2]),initial_v_m_s=float(initial[3]))]
    colors=[tuple(c) for c in plt.cm.viridis(np.linspace(.05,.8,len(cfg['diameters_um'])))]
    measured_re=dict(y_mm=y,height_above_first_measurement_mm=height,mean_u_m_s=u,mean_v_m_s=v,
                     x_mm=profile['x_mm'],N=profile['N'],native_row=profile['native_row'],native_column=profile['native_column'])
    trajectories={};checks=[];combined=[]
    gas=lambda position:np.asarray(cfg['gas_velocity_m_s'])
    for d_um in cfg['diameters_um']:
        d=d_um*1e-6;tau=float(model.relaxation_time(d,p))
        tr,status=model.upward_trajectory(d,initial,top,gas,p,cfg)
        if status=='turning_point':turning,turn_status=tr,status
        else:turning,turn_status=model.upward_trajectory(d,initial,np.inf,gas,p,cfg)
        endpoint=endpoint_record('stationary_air',d_um,initial,tr,status,turning,turn_status,y0)
        endpoints.append(endpoint)
        re=model.reynolds(d,tr[:,3:5],cfg['gas_velocity_m_s'],p)
        columns=dict(assumed_diameter_um=np.full(len(tr),d_um),t_s=tr[:,0],x_mm=tr[:,1]*1e3,
                     y_mm=tr[:,2]*1e3,height_above_first_measurement_mm=tr[:,2]*1e3-y0,
                     u_m_s=tr[:,3],v_m_s=tr[:,4],Re_in_plane=re,
                     tau_stokes_s=np.full(len(tr),tau),tau_drag_secant_s=tau/model.drag_factor(re,cfg))
        csv_name=f'trajectory_stationary_air_d{d_um:g}um.csv'
        model.write_csv(out/csv_name,columns);trajectories[d_um]=csv_name
        combined.append(columns)
        measured_re[f'Re_mean_vector_d{d_um:g}um']=model.reynolds(d,w,cfg['gas_velocity_m_s'],p)
        # Compare each endpoint variable in its own units, including actual turning point.
        dt=min(cfg['trajectory_time_step_s'],tau/30)/2
        fine,fine_status=model.upward_trajectory(d,initial,top,gas,p,cfg,dt=dt)
        fine_turn,fine_turn_status=(fine,fine_status) if fine_status=='turning_point' else model.upward_trajectory(d,initial,np.inf,gas,p,cfg,dt=dt)
        assert status==fine_status and turn_status==fine_turn_status
        delta=np.abs(tr[-1]-fine[-1]);turn_delta=np.abs(turning[-1]-fine_turn[-1])
        assert delta[0]<2e-6 and np.max(delta[1:3])<2e-7 and np.max(delta[3:])<2e-6
        assert turn_delta[0]<2e-6 and np.max(turn_delta[1:3])<2e-7
        checks.append(dict(diameter_um=d_um,dt_s=min(cfg['trajectory_time_step_s'],tau/30),
                           half_step_endpoint_time_difference_s=float(delta[0]),
                           half_step_endpoint_position_difference_m=delta[1:3].tolist(),
                           half_step_endpoint_velocity_difference_m_s=delta[3:].tolist(),
                           half_step_turning_height_difference_m=float(turn_delta[2]),max_Re=float(re.max())))
    model.write_csv(out/'reynolds_measured_velocity_assumed_air.csv',measured_re)
    model.write_csv(out/'reynolds_reference_trajectories.csv',
                    {k:np.concatenate([c[k] for c in combined]) for k in combined[0]})
    model.write_csv(out/'reference_endpoints.csv',{k:[e[k] for e in endpoints] for k in endpoints[0]})
    manifest=[]
    # Both velocity views use exactly the same native CSV curves, styles and ordinate.
    for stem,xlim in [('velocity_reference_full',(0,hmax)),
                      ('velocity_reference_near_origin',tuple(cfg['velocity_figure_near_origin_limits_mm']))]:
        fig,ax=axis(r'$h=y-y_0\ (\mathrm{mm})$',r'$v,\ \overline{v}\ (\mathrm{m\,s^{-1}})$')
        curves=[curve(ax,out,'measured_fixed_x_profile.csv','height_above_first_measurement_mm','mean_v',
                      'Ensemble mean (fixed x)',color='#d62728',linewidth=1.7),
                curve(ax,out,'ballistic_reference.csv','height_above_first_measurement_mm','v_m_s',
                      'Gravity only',color='black',linestyle='--')]
        for d_um,color in zip(cfg['diameters_um'],colors):
            curves.append(curve(ax,out,trajectories[d_um],'height_above_first_measurement_mm','v_m_s',
                                rf'$d={d_um:g}\,\mu\mathrm{{m}}$',color=color))
        # Include all signed measured values; no ordinate clipping of negative measurements.
        lower=min(0.,float(np.nanmin(v)));upper=max(float(np.nanmax(v)),initial[3]);padding=.04*max(upper-lower,.01)
        ax.set(xlim=xlim,ylim=(lower-padding,upper+padding))
        export(fig,ax,out,stem,curves,manifest,velocity=True)
    fig,ax=axis(r'$h=y-y_0\ (\mathrm{mm})$',r'$Re_d^{\mathrm{mean}}\ (-)$');curves=[]
    for d_um,color in zip(cfg['diameters_um'],colors):
        curves.append(curve(ax,out,'reynolds_measured_velocity_assumed_air.csv','height_above_first_measurement_mm',
            f'Re_mean_vector_d{d_um:g}um',rf'$d={d_um:g}\,\mu\mathrm{{m}}$',color=color))
    ax.set_xlim(0,hmax)
    export(fig,ax,out,'reynolds_measured_velocity_assumed_air',curves,manifest,category='conditional_estimate')
    fig,ax=axis(r'$h=y-y_0\ (\mathrm{mm})$',r'$Re_d^{\mathrm{ref}}\ (-)$');curves=[]
    for d_um,color in zip(cfg['diameters_um'],colors):
        curves.append(curve(ax,out,'reynolds_reference_trajectories.csv','height_above_first_measurement_mm','Re_in_plane',
            rf'$d={d_um:g}\,\mu\mathrm{{m}}$',color=color,filter_column='assumed_diameter_um',filter_value=d_um))
    ax.set_xlim(0,hmax)
    export(fig,ax,out,'reynolds_reference_trajectories',curves,manifest)
    diam=np.geomspace(min(cfg['diameters_um']),max(cfg['diameters_um']),201);tau=model.relaxation_time(diam*1e-6,p)
    model.write_csv(out/'particle_relaxation_time.csv',dict(assumed_diameter_um=diam,tau_stokes_s=tau,tau_stokes_ms=tau*1e3))
    fig,ax=axis(r'$d\ (\mu\mathrm{m})$',r'$\tau_{p,\mathrm{Stokes}}\ (\mathrm{ms})$')
    curves=[curve(ax,out,'particle_relaxation_time.csv','assumed_diameter_um','tau_stokes_ms',
                  'Stokes-limit response time',color='tab:blue')]
    ax.set(xscale='log',yscale='log');export(fig,ax,out,'particle_relaxation_time',curves,manifest)
    stop_h=np.array([e['actual_turning_height_mm'] for e in endpoints[1:]])
    end_h=np.array([e['plotted_endpoint_height_mm'] for e in endpoints[1:]])
    sensitivity=dict(assumed_diameter_um=cfg['diameters_um'],actual_stationary_air_turning_height_mm=stop_h,
        plotted_stationary_air_endpoint_height_mm=end_h,plotted_endpoint_reason=[e['plotted_endpoint_reason'] for e in endpoints[1:]],
        tau_stokes_s=[float(model.relaxation_time(d*1e-6,p)) for d in cfg['diameters_um']],
        gravity_only_actual_turning_height_mm=np.full(len(stop_h),h_stop),measurement_height_extent_mm=np.full(len(stop_h),hmax))
    model.write_csv(out/'diameter_sensitivity.csv',sensitivity)
    fig,ax=axis(r'$d\ (\mu\mathrm{m})$',r'$h_{\mathrm{turn}}\ (\mathrm{mm})$')
    curves=[curve(ax,out,'diameter_sensitivity.csv','assumed_diameter_um','actual_stationary_air_turning_height_mm',
                  'Stationary-air turning height',color='tab:blue'),
            curve(ax,out,'diameter_sensitivity.csv','assumed_diameter_um','gravity_only_actual_turning_height_mm',
                  'Gravity-only turning height',color='black',linestyle='--'),
            curve(ax,out,'diameter_sensitivity.csv','assumed_diameter_um','measurement_height_extent_mm',
                  'Measurement boundary',color='0.5',linestyle=':')]
    ax.set_yscale('log');export(fig,ax,out,'diameter_sensitivity',curves,manifest,category='supplementary')
    jet_input=model.ROOT/cfg['znmf_profile_csv'];z=read_csv(jet_input)
    # String status is checked separately; genfromtxt numeric reader above preserves numerical columns.
    structured=np.genfromtxt(jet_input,names=True,delimiter=',',dtype=None,encoding='utf-8')
    keep=(structured['status']=='measured')&np.isfinite(z['Uc_m_per_s'])&(z['Uc_m_per_s']>0)
    jet_x=z['s_mm'][keep];uc=z['Uc_m_per_s'][keep]
    jet_columns=dict(x_ZNMF_mm=jet_x,Uc_m_s=uc,tau_flow_D0_over_Uc_s=cfg['znmf_nozzle_diameter_m']/uc)
    for d_um in cfg['diameters_um']:
        jet_columns[f'St_D_d{d_um:g}um']=model.relaxation_time(d_um*1e-6,p)*uc/cfg['znmf_nozzle_diameter_m']
    jet_csv='Jet_Interaction_Sensitivity/znmf_local_nozzle_stokes.csv';model.write_csv(out/jet_csv,jet_columns)
    fig,ax=axis(r'$x_{\mathrm{ZNMF}}\ (\mathrm{mm})$',r'$St_D(x)\ (-)$');curves=[]
    for d_um,color in zip(cfg['diameters_um'],colors):
        curves.append(curve(ax,out,jet_csv,'x_ZNMF_mm',f'St_D_d{d_um:g}um',rf'$d={d_um:g}\,\mu\mathrm{{m}}$',color=color))
    ax.set_yscale('log');export(fig,ax,out,FIGURES[-1],curves,manifest,category='hypothetical_jet_interaction')
    summary=dict(initial_y_mm=y0,initial_u_m_s=float(initial[2]),initial_v_m_s=float(initial[3]),
                 initial_native_row=int(profile['native_row'][first]),initial_native_column=int(profile['native_column'][first]),
                 fixed_x_mm=maximum['x_mm'],measurement_height_extent_mm=hmax,maximum=maximum,endpoints=endpoints)
    params=dict(configuration=cfg,properties=p,calibration=dict(um_per_pixel=5.35,pulse_separation_s=25e-6,
        velocity_conversion_m_s_per_pixel=.214,action='Saved final SI fields used unchanged; no repeated calibration'),
        initial_condition=summary,jet_input=dict(path=str(jet_input),sha256=digest(jet_input),
        registration='Separate ZNMF coordinate retained; no spatial registration to PZT',
        input_subset='Existing separately measured comparison table; only finite positive measured Uc rows'),
        reviewed_commit='1522f54161c7',subsequent_change_check='No later commits; active source/config matched reviewed checkout before revision')
    (out/'parameters_used.json').write_text(json.dumps(params,indent=2),encoding='utf-8')
    (jet/'parameters_used.json').write_text(json.dumps(dict(definition='St_D(x)=tau_p,Stokes*Uc(x)/D0',
        D0_m=cfg['znmf_nozzle_diameter_m'],assumed_diameters_um=cfg['diameters_um'],
        tau_definition='rho_water*d^2/(18*mu_air), Stokes-limit baseline',
        provenance=params['jet_input'],interpretation='Hypothetical local nozzle-length response comparison; neither PZT-plume nor nozzle-exit Stokes number'),indent=2),encoding='utf-8')
    for name,sha in protected.items():assert digest(model.FINAL/name)==sha,name
    checks_record=dict(physics=model.validate_physics(p,cfg),calibration_once=True,SI_integration=True,
                       step_convergence=checks,experimental_files_unchanged=len(protected))
    (validation_dir/'numerical_checks.json').write_text(json.dumps(checks_record,indent=2))
    (validation_dir/'experimental_input_hashes.json').write_text(json.dumps(protected,indent=2))
    (out/'figure_curve_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    from .theory_report import write_report
    write_report(out,params,manifest)
    (out/'regenerate.py').write_text('''from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[4]))
from PZT_Atomiser_Analysis.analysis.theory_comparison import run
run()
''')
    return dict(initial_y_mm=y0,initial_v_m_s=float(initial[3]),fixed_x_mm=maximum['x_mm'],
                gravity_only_actual_turning_height_mm=h_stop,measurement_height_extent_mm=hmax,
                stationary_air_turning_heights_mm=stop_h.tolist(),figure_count=len(manifest),
                experimental_files_unchanged=len(protected),assumed_diameters_um=cfg['diameters_um'])


def notebook_cells():
    import nbformat as nb
    text='''## Fixed-line ensemble means and stationary-air trajectory references

The measured profile is an ensemble mean along a fixed vertical line, not a
trajectory or tracked centreline. References use assumed 5, 10, 20, 50 and
100 micrometre water droplets in stationary air. Height h=y-y0 starts at the
first finite paired measurement used for initialisation, not the atomiser.
The full and near-origin views use identical curves. Mean-vector conditional
Reynolds estimates and trajectory Reynolds numbers are shown separately.
The hypothetical ZNMF Stokes comparison is in Jet_Interaction_Sensitivity;
it supplies no measured gas velocity for the PZT experiment.
'''
    code='''from PZT_Atomiser_Analysis.analysis.theory_comparison import run as run_theory
from PZT_Atomiser_Analysis.analysis.theory_package import FIGURES
theory=run_theory()
print(json.dumps(theory,indent=2))
for stem in FIGURES:
    print(stem)
    display(Image(filename=str(output/'Theory_Comparison'/(stem+'.png')),width=650))
'''
    cells=[nb.v4.new_markdown_cell(text),nb.v4.new_code_cell(code)]
    for c in cells:c.metadata['theory_comparison']=True
    return cells
