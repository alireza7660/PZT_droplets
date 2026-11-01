"""Spherical-droplet references from final SI data; no experimental files rewritten.

Run from workspace root: python -m PZT_Atomiser_Analysis.analysis.theory_comparison
All integrations use SI units; plotted heights alone are converted to mm.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .plotting import plt
from .vertical_profile import extract_vertical_profile

ROOT=Path(__file__).resolve().parents[2]
PROJECT=ROOT/'PZT_Atomiser_Analysis'
FINAL=PROJECT/'Final_Results/Batch4_Cam1'
SOURCES={
 'water_density': 'https://www.nist.gov/system/files/documents/2017/05/09/NISTIR_7383_2013-04-24_2015-07-13Rev.pdf',
 'air_density': 'https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=104&refer=ap',
 'air_viscosity': 'https://ntrs.nasa.gov/api/citations/20250005707/downloads/main.pdf?attachment=true',
 'drag_and_response': 'https://doc.comsol.com/6.4/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.36.html'}


def fluid_properties(config):
    """20 C water density (NIST GLP10 polynomial), dry-air reference and Sutherland mu."""
    t=float(config['temperature_C'])
    if t!=20: raise ValueError('Re-source density properties before using another temperature')
    # NIST GLP10 air-free water, Tanaka equation at atmospheric pressure.
    water=999.974950*(1-((t-3.983035)**2*(t+301.797))/(522528.9*(t+69.34881)))
    # Rounded dry-air reference used as an explicitly configured property approximation.
    air=1.20479*(config['pressure_Pa']/101325.)
    if config['relative_humidity'] not in (None,0): raise ValueError('Humid-air properties require a supplied property model')
    mu=1.716e-5*((t+273.15)/273.15)**1.5*(273.15+110.4)/(t+273.15+110.4)
    return dict(rho_water_kg_m3=water,rho_air_kg_m3=air,mu_air_Pa_s=mu,g_m_s2=9.80665,
                air_density_note='NIST dry-air reference 1.20479 kg/m3; approximate 20 C, pressure scaled; not measured density',
                mu_air_note='Sutherland: mu0=1.716e-5 Pa s, T0=273.15 K, S=110.4 K',sources=SOURCES)


def relaxation_time(d_m,properties):
    d=np.asarray(d_m,float)
    if np.any(d<=0): raise ValueError('Diameter must be positive, in metres')
    return properties['rho_water_kg_m3']*d*d/(18*properties['mu_air_Pa_s'])


def reynolds(d_m,droplet_velocity,gas_velocity,properties):
    """Two-component in-plane slip Re; no unavailable out-of-plane component invented."""
    return properties['rho_air_kg_m3']*np.asarray(d_m)*np.linalg.norm(
        np.asarray(gas_velocity)-np.asarray(droplet_velocity),axis=-1)/properties['mu_air_Pa_s']


def drag_factor(re,config):
    """C_D Re/24, with exact Stokes branch and finite zero-slip limit."""
    re=np.asarray(re,float)
    if np.any(re<0) or np.any(re>=config['schiller_naumann_re_limit_exclusive']):
        raise ValueError('Outside configured spherical drag range 0<=Re<800')
    return np.where(re<=config['stokes_drag_re_threshold'],1.,1.+.15*re**.687)


def acceleration(velocity,d_m,gas_velocity,properties,config,drag=True,buoyancy=True):
    slip=np.asarray(gas_velocity,float)-np.asarray(velocity,float)
    a=np.zeros_like(slip)
    a[...,1]=-properties['g_m_s2']*(1-properties['rho_air_kg_m3']/properties['rho_water_kg_m3'] if buoyancy else 1)
    if drag:
        re=reynolds(d_m,velocity,gas_velocity,properties)
        a=a+slip*drag_factor(re,config)/relaxation_time(d_m,properties)
    return a


def upward_trajectory(d_m,initial,top_m,gas,properties,config,dt=None):
    """RK4 in time; bracket/bisect first turning or top-domain event.

    initial=(x_m,y_m,u_m_s,v_m_s). gas(position_m) returns (u_air,v_air).
    Returns ascending branch only. No empirical-profile interpolation is used.
    """
    state=np.asarray(initial,float).copy()
    if state[3]<=0: raise ValueError('Upward reference needs a positive measured initial v')
    step=min(dt or config['trajectory_time_step_s'],float(relaxation_time(d_m,properties))/30)
    def rhs(s): return np.r_[s[2:],acceleration(s[2:],d_m,gas(s[:2]),properties,config)]
    def advance(s,h):
        a=rhs(s);b=rhs(s+h*a/2);c=rhs(s+h*b/2);d=rhs(s+h*c)
        return s+h*(a+2*b+2*c+d)/6
    rows=[np.r_[0.,state]];time=0.;status='time_limit'
    while time<config['trajectory_max_time_s']:
        h=min(step,config['trajectory_max_time_s']-time);new=advance(state,h)
        events=[]
        for name,index,target in [('turning_point',3,0.),('domain_top',1,top_m)]:
            crossed=(new[index]<=target if index==3 else new[index]>=target)
            if crossed:
                lo,hi=0.,h
                for _ in range(42):
                    mid=(lo+hi)/2;s=advance(state,mid)
                    if (s[index]<=target if index==3 else s[index]>=target): hi=mid
                    else:lo=mid
                events.append((hi,name,index,target))
        if events:
            h,status,index,target=min(events);new=advance(state,h);new[index]=target
        time+=h;rows.append(np.r_[time,new]);state=new
        if events:break
    result=np.asarray(rows)
    if np.any(np.diff(result[:,2])<-1e-12) or np.any(result[:,4]<0):raise AssertionError('Descending branch')
    return result,status


def write_csv(path,columns):
    keys=list(columns);arrays=[np.asarray(columns[k]) for k in keys]
    if len({len(a) for a in arrays})!=1:raise ValueError('Unaligned CSV columns')
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(keys);writer.writerows(zip(*arrays))


def figure(out,name,xlabel,ylabel):
    fig,ax=plt.subplots(figsize=(6.4,4.5),layout='constrained');ax.set(xlabel=xlabel,ylabel=ylabel)
    return fig,ax


def export(fig,ax,out,name):
    if name=='velocity_theory_comparison':
        ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2)
    else:ax.legend(frameon=False,fontsize=8)
    for ext in ('pdf','svg','png'):fig.savefig(out/(name+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)


def load_gas(config,profile):
    """Optional simultaneous gas data must be native-y aligned; never auto-register."""
    path=config.get('measured_gas_profile_csv')
    if not path:return None
    if not config.get('measured_gas_profile_registration'):raise ValueError('Document simultaneous gas registration')
    z=np.genfromtxt(ROOT/path,names=True,delimiter=',')
    np.testing.assert_allclose(z['y_mm'],profile['y_mm'],atol=1e-10,rtol=0)
    a=np.c_[z['u_air_m_s'],z['v_air_m_s']]
    if not np.isfinite(a).all():raise ValueError('Registered gas profile contains gaps')
    return a


def notebook_cells():
    import nbformat as nb
    text='''## Droplet theory and dimensionless sensitivity references

The final calibrated droplet data are used unchanged. Diameter and simultaneous gas velocity are not measured inputs in the available files. Diameter cases and still-air references are assumptions. The fixed-x ensemble mean is not a trajectory. ZNMF is a separate gas experiment used only for an explicitly hypothetical jet-interaction Stokes comparison. Height is h=y-y0 relative to the first finite paired measurement, not the atomiser surface. See Theory_Comparison/README.md and methods_results.tex for definitions and limitations.
'''
    code='''from PZT_Atomiser_Analysis.analysis.theory_comparison import run as run_theory
theory = run_theory()
print(json.dumps(theory, indent=2))
for name in ['velocity_theory_comparison','reynolds_vs_height','particle_relaxation_time','diameter_sensitivity','znmf_local_nozzle_stokes']:
    display(Image(filename=str(output/'Theory_Comparison'/(name+'.png')),width=650))
'''
    cells=[nb.v4.new_markdown_cell(text),nb.v4.new_code_cell(code)]
    for cell in cells:cell.metadata['theory_comparison']=True
    return cells


def run(config_path=None):
    cfg_path=Path(config_path or PROJECT/'config/theory_comparison.json')
    cfg=json.loads(cfg_path.read_text());props=fluid_properties(cfg)
    if np.any(np.asarray(cfg['gas_velocity_m_s'])!=0):
        raise ValueError('This reference is labelled quiescent air; use a separately labelled scenario for nonzero gas')
    out=FINAL/'Theory_Comparison';out.mkdir(exist_ok=True)
    # Snapshot all existing experimental numerical and figure files, excluding new model directory.
    protected=[p for p in FINAL.rglob('*') if p.is_file() and out not in p.parents and p.suffix.lower() in ('.npz','.csv','.pdf','.svg','.png')]
    hashes={str(p.relative_to(FINAL)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    active=json.loads((PROJECT/'config/analysis.json').read_text())
    np.testing.assert_allclose([active['coordinate_scale'],active['velocity_scale'],active['confirmed_pulse_separation_seconds']],[.00535,.214,25e-6],rtol=0,atol=1e-15)
    profile,maximum=extract_vertical_profile(FINAL/'data/statistics.npz')
    with np.load(FINAL/'data/statistics.npz') as z:
        u=z['mean_u'][profile['native_row'],profile['native_column']]
        keys=z.files
    y=profile['y_mm'];v=profile['mean_v'];velocity=np.c_[u,v]
    good=np.isfinite(velocity).all(axis=1)&(profile['N']>=cfg['minimum_initial_sample_count'])
    first=np.flatnonzero(good)[0]
    initial=np.array([profile['x_mm'][first]*1e-3,y[first]*1e-3,u[first],v[first]])
    y0=float(y[first]);height=y-y0
    props_record=dict(properties=props,configuration=cfg,source_configuration=str(cfg_path),
        calibration=dict(um_per_pixel=5.35,dt_s=25e-6,velocity_factor_m_s_per_pixel=.214,action='Read final SI velocities unchanged; no calibration reapplied'),
        initial_condition=dict(native_row=int(profile['native_row'][first]),native_column=int(profile['native_column'][first]),state_SI=initial.tolist(),height_reference_y0_mm=y0,comparison_height_definition='h_mm=y_mm-y0_mm; not distance from atomiser surface',reliability='Finite paired means, N>=configured minimum only; no physical validation claim'),
        maximum=maximum)
    import h5py
    with h5py.File(ROOT/'RAW_data/Batch4_Cam1.h5','r') as h5:
        fields=h5['PIV FIELDS'];first_key=next(iter(fields))
        inventory=dict(hdf5_path='RAW_data/Batch4_Cam1.h5',root_groups=list(h5),root_attribute_names=list(h5.attrs),
            field_count=len(fields),representative_field=first_key,field_datasets=list(fields[first_key]),
            statistics_datasets=list(h5['STATISTICS']),saved_statistic_arrays=keys,
            instantaneous_final_arrays=['frame_1.npz: u,v','frame_879.npz: u,v','frame_1757.npz: u,v'],
            velocity_method='PIV image-correlation export, not identified droplet tracks',
            diameter_evidence='No validated sizing result located in active project or RAW_data; sensitivity only',
            simultaneous_gas='Not found; separate ZNMF experiment exists',
            spatial_origin=active['origin_definition'])
    (out/'input_inventory.json').write_text(json.dumps(inventory,indent=2))
    measured=load_gas(cfg,profile)
    write_csv(out/'measured_fixed_x_profile.csv',dict(**profile,mean_u=u,height_above_first_measurement_mm=height))
    fig,ax=figure(out,'velocity_theory_comparison',r'$h=y-y_0\ (\mathrm{mm})$',r'$v,\ \overline{v}\ (\mathrm{m\,s^{-1}})$')
    ax.plot(height,v,color='#d62728',lw=1.5,label=rf'Measured fixed $x={maximum["x_mm"]:.5f}$ mm')
    # Gravity-only ballistic reference: no aerodynamic force, no buoyancy.
    g=props['g_m_s2'];end=min(y.max()*1e-3,initial[1]+initial[3]**2/(2*g))
    by=np.linspace(initial[1],end,501);bv=np.sqrt(np.maximum(0,initial[3]**2-2*g*(by-initial[1])))
    ax.plot(by*1e3-y0,bv,'k--',lw=1,label='Gravity-only ballistic reference')
    write_csv(out/'ballistic_reference.csv',dict(y_mm=by*1e3,height_above_first_measurement_mm=by*1e3-y0,v_m_s=bv))
    colors=plt.cm.viridis(np.linspace(.05,.8,len(cfg['diameters_um'])))
    cases=[];re_cols=dict(y_mm=y,height_above_first_measurement_mm=height,mean_u_m_s=u,mean_v_m_s=v)
    endpoints=[]
    for d_um,color in zip(cfg['diameters_um'],colors):
        d=d_um*1e-6
        gas=lambda position:np.asarray(cfg['gas_velocity_m_s'])
        trajectory,status=upward_trajectory(d,initial,y.max()*1e-3,gas,props,cfg)
        re=reynolds(d,trajectory[:,3:5],cfg['gas_velocity_m_s'],props)
        tau=float(relaxation_time(d,props));factor=drag_factor(re,cfg)
        write_csv(out/f'trajectory_quiescent_d{d_um:g}um.csv',dict(t_s=trajectory[:,0],x_mm=trajectory[:,1]*1e3,y_mm=trajectory[:,2]*1e3,height_above_first_measurement_mm=trajectory[:,2]*1e3-y0,u_m_s=trajectory[:,3],v_m_s=trajectory[:,4],Re_in_plane=re,tau_stokes_s=np.full(len(re),tau),tau_drag_secant_s=tau/factor))
        ax.plot(trajectory[:,2]*1e3-y0,trajectory[:,4],color=color,label=rf'Assumed still air, $d={d_um:g}\,\mu$m')
        re_cols[f'Re_still_air_d{d_um:g}um']=reynolds(d,velocity,cfg['gas_velocity_m_s'],props)
        if measured is not None:re_cols[f'Re_measured_gas_d{d_um:g}um']=reynolds(d,velocity,measured,props)
        # Integration convergence: repeat with half step, compare event height/time.
        fine,fine_status=upward_trajectory(d,initial,y.max()*1e-3,gas,props,cfg,dt=min(cfg['trajectory_time_step_s'],tau/30)/2)
        error=float(np.max(np.abs(trajectory[-1]-fine[-1])))
        if error>2e-6 or status!=fine_status:raise AssertionError('Trajectory step convergence failed')
        cases.append(dict(diameter_um=d_um,tau_stokes_s=tau,termination=status,end_y_mm=float(trajectory[-1,2]*1e3),end_t_s=float(trajectory[-1,0]),max_Re=float(re.max()),half_step_max_state_difference=error))
        endpoints.append(trajectory[-1,2]*1e3)
    export(fig,ax,out,'velocity_theory_comparison')
    write_csv(out/'reynolds_vs_height.csv',re_cols)
    fig,ax=figure(out,'reynolds_vs_height',r'$h=y-y_0\ (\mathrm{mm})$',r'$Re_d\ (-)$')
    for d_um,color in zip(cfg['diameters_um'],colors):
        ax.plot(height,re_cols[f'Re_still_air_d{d_um:g}um'],color=color,label=rf'Still air, assumed $d={d_um:g}\,\mu$m')
        if measured is not None:ax.plot(height,re_cols[f'Re_measured_gas_d{d_um:g}um'],'--',color=color,label=f'Measured gas, assumed d={d_um:g} um')
    export(fig,ax,out,'reynolds_vs_height')
    diam=np.geomspace(min(cfg['diameters_um']),max(cfg['diameters_um']),201)
    tau=relaxation_time(diam*1e-6,props);st=dict(diameter_um=diam,tau_stokes_s=tau)
    fig,ax=figure(out,'particle_relaxation_time',r'$d\ (\mu\mathrm{m})$',r'$\tau_p\ (\mathrm{ms})$')
    ax.loglog(diam,tau*1e3,label='Stokes-limit spherical-droplet response time')
    export(fig,ax,out,'particle_relaxation_time');write_csv(out/'particle_relaxation_time.csv',st)
    # Preserve the superseded arbitrary-timescale references outside active outputs.
    for ext in ('pdf','svg','png','csv'):
        old=out/('stokes_vs_diameter.'+ext)
        if old.exists():
            archive=out/'archive_assumed_flow_timescale';archive.mkdir(exist_ok=True)
            old.replace(archive/old.name)
    fig,ax=figure(out,'diameter_sensitivity',r'$d\ (\mu\mathrm{m})$',r'$y_{\mathrm{end}}-y_0\ (\mathrm{mm})$')
    rise=np.asarray(endpoints)-initial[1]*1e3
    ax.plot(cfg['diameters_um'],rise,'o-',color='tab:blue',label='Assumed still-air ascending reference')
    ax.axhline((end-initial[1])*1e3,color='k',ls='--',label='Ballistic rise within domain')
    export(fig,ax,out,'diameter_sensitivity')
    write_csv(out/'diameter_sensitivity.csv',dict(diameter_um=cfg['diameters_um'],ascending_rise_within_domain_mm=rise,ballistic_rise_within_domain_mm=np.full(len(rise),(end-initial[1])*1e3),tau_stokes_s=[c['tau_stokes_s'] for c in cases],max_Re=[c['max_Re'] for c in cases],termination=[c['termination'] for c in cases]))
    znmf=ROOT/cfg['znmf_profile_csv']
    z=np.genfromtxt(znmf,names=True,delimiter=',',dtype=None,encoding='utf-8')
    status=z['status']=='measured';positive=status&np.isfinite(z['Uc_m_per_s'])&(z['Uc_m_per_s']>0)
    coord=z['s_mm'][positive];uc=z['Uc_m_per_s'][positive]
    zcols=dict(znmf_own_axial_coordinate_mm=coord,Uc_m_s=uc)
    fig,ax=figure(out,'znmf_local_nozzle_stokes',r'$x_{\mathrm{ZNMF}}\ (\mathrm{mm})$',r'$St_D(x)=\tau_p U_c(x)/D_0\ (-)$')
    for d_um,color in zip(cfg['diameters_um'],colors):
        number=relaxation_time(d_um*1e-6,props)*uc/cfg['znmf_nozzle_diameter_m']
        zcols[f'St_D_d{d_um:g}um']=number
        ax.plot(coord,number,color=color,label=rf'Hypothetical jet, $d={d_um:g}\,\mu$m')
    export(fig,ax,out,'znmf_local_nozzle_stokes');write_csv(out/'znmf_local_nozzle_stokes.csv',zcols)
    write_csv(out/'diameter_case_summary.csv',dict(
        assumed_diameter_um=cfg['diameters_um'],tau_stokes_s=[c['tau_stokes_s'] for c in cases],
        min_Re_measured_mean_still_air=[float(np.nanmin(re_cols[f'Re_still_air_d{d:g}um'])) for d in cfg['diameters_um']],
        max_Re_measured_mean_still_air=[float(np.nanmax(re_cols[f'Re_still_air_d{d:g}um'])) for d in cfg['diameters_um']],
        min_hypothetical_jet_St_D=[float(np.min(zcols[f'St_D_d{d:g}um'])) for d in cfg['diameters_um']],
        max_hypothetical_jet_St_D=[float(np.max(zcols[f'St_D_d{d:g}um'])) for d in cfg['diameters_um']],
        quiescent_reference_rise_mm=rise))
    props_record['znmf_input']=dict(path=str(znmf),sha256=hashlib.sha256(znmf.read_bytes()).hexdigest(),registration='None; separate experiment coordinate; no combined-flow calculation')
    (out/'parameters_used.json').write_text(json.dumps(props_record,indent=2),encoding='utf-8')
    checks=validate_physics(props,cfg)
    unchanged=all(hashlib.sha256((FINAL/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    if not unchanged:raise AssertionError('Existing experimental output changed')
    summary=dict(cases=cases,maximum=maximum,initial_y_mm=float(initial[1]*1e3),comparison_height_range_mm=[float(height.min()),float(height.max())],initial_mean_v_m_s=float(initial[3]),
                 fixed_x_mm=maximum['x_mm'],protected_files_unchanged=len(hashes),physics_checks=checks,
                 measured_gas_available=measured is not None,diameters_measured=False)
    (out/'validation.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (out/'experimental_input_hashes.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
    (out/'regenerate.py').write_text('''"""Run the reusable comparison from this results folder or any working directory."""
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(root))
from PZT_Atomiser_Analysis.analysis.theory_comparison import run
print(run())
''',encoding='utf-8')
    write_methods(out,cfg,props,summary)
    return summary


def validate_physics(p,c):
    d=20e-6;zero=np.zeros(2);g=p['g_m_s2']*(1-p['rho_air_kg_m3']/p['rho_water_kg_m3'])
    np.testing.assert_allclose(acceleration(zero,d,zero,p,c),[0,-g])
    up=acceleration([0,1],d,zero,p,c);down=acceleration([0,-1],d,zero,p,c)
    assert up[1]<-g and down[1]>-g
    tiny=np.array([1e-9,2e-9]);a=acceleration(tiny,d,zero,p,c)-np.array([0,-g])
    np.testing.assert_allclose(a,-tiny/relaxation_time(d,p),rtol=1e-6,atol=1e-14)
    np.testing.assert_allclose(relaxation_time(2*d,p),4*relaxation_time(d,p))
    np.testing.assert_allclose(reynolds(d,[3,4],[0,0],p),p['rho_air_kg_m3']*d*5/p['mu_air_Pa_s'])
    try:drag_factor(800,c)
    except ValueError:pass
    else:raise AssertionError('Drag validity not enforced')
    tiny_d=5e-6;v0=.001;tau=float(relaxation_time(tiny_d,p))
    trajectory,status=upward_trajectory(tiny_d,[0,0,0,v0],1.,lambda q:zero,p,c)
    analytic_time=tau*np.log1p(v0/(g*tau))
    analytic_height=tau*v0-g*tau**2*np.log1p(v0/(g*tau))
    assert status=='turning_point'
    np.testing.assert_allclose(trajectory[-1,[0,2]],[analytic_time,analytic_height],rtol=1e-6,atol=1e-14)
    return dict(zero_slip_finite=True,gravity_downward=True,drag_opposes_slip=True,
                low_Re_Stokes_limit=True,analytic_Stokes_turning_time_height=True,tau_scales_as_d_squared=True,vector_norm_3_4_5=True,Re_validity_enforced=True)


def write_methods(out,cfg,p,summary):
    text=f'''# Droplet theory comparison: conditional references

Final PZT data: 190 kHz, 9 Vpp; 5.35 micrometres/pixel, 25 microseconds, 0.214 m/s per pixel displacement. Saved SI fields are read without applying calibration again. x/u are horizontal/rightward, y/v upward. All signed experimental velocities remain unchanged.

## Available evidence and missing inputs

`../data/statistics.npz` provides ensemble mean u/v, RMS, covariance and N; `../data/frame_1.npz`, frame_879 and frame_1757 provide representative instantaneous vector fields. RAW_data/Batch4_Cam1.h5 contains correlation-export fields us/vs and flags, not droplet identities or trajectories. analysis/calibrated.py accumulates finite paired samples; analysis/final_analysis.py adopts the final convention. No validated sizing result or simultaneous gas field was found in the active project or RAW_data. Image-spot width, the water-flow-rate image and nozzle dimensions are not diameter evidence. The exported (0,0) reference is not a located atomiser surface. The adopted component/displacement interpretation remains documented in the experimental configuration; this comparison does not independently verify it.

The fixed-x profile is column 42 through the global mean-v maximum, not a tracked centreline. A declining value there may partly reflect lateral plume displacement. Neither a fixed-position profile nor a height-wise maximum of ensemble means is an individual droplet trajectory. No height-wise maximum is substituted here.

## Properties and model

Air temperature is recorded as 20 C. Pressure 101325 Pa, water also at 20 C, and constant spherical diameter are assumptions. Humidity is unavailable; dry-air properties are an approximation, not a recorded zero-humidity condition. Evaporation is neglected. rho_air={p['rho_air_kg_m3']:.8g} kg/m3; rho_water={p['rho_water_kg_m3']:.9g} kg/m3; mu_air={p['mu_air_Pa_s']:.9g} Pa s; g=9.80665 m/s2. Water density uses NIST GLP10's Tanaka polynomial; air uses the approximate NIST dry-air reference; Sutherland viscosity uses NASA's tabulated constants. No liquid viscosity enters Re or relaxation time.

Vector equations use dq/dt=w_d and dw_d/dt=(0,-g(1-rho_air/rho_water))+(w_air-w_d) f(Re)/tau_p, where Re=rho_air*d*norm(w_air-w_d)/mu_air, tau_p=rho_water*d^2/(18*mu_air). f=1 for Re<=0.01 (creeping-flow approximation); otherwise f=1+0.15 Re^0.687, equivalent to C_D=24*f/Re. This implementation enforces Re<800 following the cited COMSOL range. Using slip*f/tau avoids division by Re at zero slip. The rigid-sphere approximation omits internal circulation, deformation, evaporation, added mass, history forces, collisions and two-way coupling. These may limit quantitative applicability.

Both in-plane components enter drag and Re; the unmeasured out-of-plane component is omitted, so Re is an in-plane-slip estimate. Re of a mean vector is not the ensemble mean of instantaneous Re. No instantaneous sample population is relabelled as a measured diameter-resolved quantity. Assumed diameters {cfg['diameters_um']} micrometres are sensitivity cases, not an inferred size distribution.

The initial state is the lowest finite paired mean-u/v sample with N>={cfg['minimum_initial_sample_count']}, y0={summary['initial_y_mm']:.8g} mm, v0={summary['initial_mean_v_m_s']:.9g} m/s at fixed x={summary['fixed_x_mm']:.8g} mm. This is a sample-availability criterion, not independent physical validation. It is not an outlet launch condition. All PZT comparison heights use h=y-y0; absolute image-reference y values are also retained in CSVs. RK4 integrates individual monodisperse reference trajectories in assumed stationary air, including horizontal motion, and stops at the first v=0 turning point or measured upper boundary. Each reference starts at the measured x but may drift horizontally: it is not a prediction along the fixed-x measurement column. Gravity-only ballistic motion excludes both drag and buoyancy. No descending branch is joined to the ascending curve; no terminal settling prediction is compared to upward plume velocity.

## Response question and Stokes number

The particle relaxation time is defined as the Stokes-limit baseline tau_p=rho_water*d^2/(18*mu_air). It is reported for every assumed diameter, even when finite-Re drag is needed; it is not represented as an exact finite-Re relaxation time. Along trajectories tau_drag_secant=tau_p/f(Re) is also exported: this is an instantaneous drag/slip coefficient timescale, not a linear perturbation response time. The longitudinal linearized drag response would involve the derivative of slip*f. No arbitrary gas disturbance duration or atomiser-frequency timescale is assigned to the active Stokes-number results.

For an explicitly hypothetical jet-interaction comparison, the separate ZNMF table supplies its measured local Uc(s) for St_D(s)=tau_p*Uc(s)/0.002. This is a local nozzle-length-based comparison, not nozzle-exit Stokes number, and not a measured PZT Stokes number. Its own s origin is retained without registering it to PZT y. No combined-flow or entrainment claim is made and no ZNMF velocity is substituted into PZT slip. The 190 kHz atomiser frequency is not used as a gas-flow frequency. No frequency-based Stokes number is computed.

## Results and limits

The measured maximum remains 1.686548142 m/s at x=14.54665 mm, y=1.284 mm. The quiescent-air model curves describe conditional deceleration and rise distances for the selected diameters, not fitted or validated experimental trajectories. Diameter, coflow and spatial sampling can all alter comparison; disagreement cannot uniquely identify any one cause. Without measured diameter and simultaneous gas velocity there is no uniquely measured Re_d or PZT Stokes number. The separate ZNMF profile exceeding some droplet velocities does not establish successful entrainment.

## Files and regeneration

Every plotted line has CSV values: measured_fixed_x_profile.csv, ballistic_reference.csv, trajectory_quiescent_d*um.csv, reynolds_vs_height.csv, particle_relaxation_time.csv, diameter_sensitivity.csv and znmf_local_nozzle_stokes.csv. parameters_used.json records properties, assumptions, sources, initial state and calibration provenance. validation.json records signs, units/scale checks, low-Re limit, numerical step convergence and unchanged experimental hashes. Configuration: ../../../../config/theory_comparison.json. Reusable source: ../../../../analysis/theory_comparison.py.

From workspace root: `.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.theory_comparison`.

Optional simultaneous gas data must use identical native y nodes and documented registration in the configuration. Its two components can produce measured-gas-slip Re sensitivity cases. A spatial trajectory comparison in a nonuniform measured gas field additionally needs a spatially resolved gas interpolant valid along the droplet path; the fixed-x gas column alone is not used to invent one.

## Figure captions

velocity_theory_comparison: Signed fixed-x mean axial droplet velocity versus h=y-y0 relative to the first finite paired measurement, gravity-only ballistic reference and assumed quiescent-air spherical-droplet trajectories for configured diameters. Model branches terminate at their turning points or the upper domain; the ensemble profile is not a trajectory.

reynolds_vs_height: In-plane-slip Reynolds estimates versus h=y-y0 from the measured mean droplet vector under assumed stationary gas and configured diameters. No measured diameter or gas velocity is implied.

particle_relaxation_time: Stokes-limit particle relaxation time versus assumed constant diameter; finite-Re drag requires the separately documented correction.

diameter_sensitivity: Ascending rise within the measured domain for assumed still-air droplets; endpoints are turning points unless marked domain_top in the CSV. The ballistic reference uses the same measured initial axial velocity.

znmf_local_nozzle_stokes: Hypothetical jet-interaction local nozzle-length Stokes comparison using separately measured ZNMF Uc(s), D0=2 mm, and assumed droplet diameters. No PZT gas registration or coflow is implied.

## Sources
'''
    for label,url in SOURCES.items():text+=f'\n- {label}: {url}\n'
    text+='\n## Computed sensitivity summary\n\n| Assumed d (micrometres) | Stokes-limit tau (ms) | Rise from y0 (mm) | Maximum trajectory in-plane Re |\n|---|---|---|---|\n'
    for case in summary['cases']:
        text+=f"| {case['diameter_um']} | {case['tau_stokes_s']*1e3:.6g} | {case['end_y_mm']-summary['initial_y_mm']:.6g} | {case['max_Re']:.6g} |\n"
    text=text.replace('../../../../config/','../../../config/').replace('../../../../analysis/','../../../analysis/')
    text=text.replace("NASA's tabulated constants","constants reported in NASA's wake-validation paper")
    peak=summary['maximum']
    text=text.replace('column 42',f"column {peak['native_column']}")
    text=text.replace('1.686548142 m/s at x=14.54665 mm, y=1.284 mm',f"{peak['mean_v_m_s']:.10g} m/s at x={peak['x_mm']:.8g} mm, y={peak['y_mm']:.8g} mm")
    (out/'README.md').write_text(text,encoding='utf-8')
    tex=r'''\section{Conditional droplet dynamics comparison}
Final calibrated PZT correlation-derived velocities are used unchanged at
$5.35\,\mu\mathrm{m/pixel}$ and $\Delta t=25\,\mu\mathrm{s}$.
The atomiser operates at 190 kHz and 9 Vpp. Coordinates $x,u$ are horizontal
and $y,v$ are upward. The image origin is not an established atomiser surface.
For a constant-diameter spherical water droplet,
\begin{align}
\dot{\boldsymbol q}&=\boldsymbol w_d,\\
\dot{\boldsymbol w}_d&=-g(1-\rho_a/\rho_w)\boldsymbol e_y+
\frac{f(Re_d)}{\tau_p}(\boldsymbol w_a-\boldsymbol w_d),\\
Re_d&=\frac{\rho_a d\lVert\boldsymbol w_a-\boldsymbol w_d\rVert}{\mu_a},
\qquad \tau_p=\frac{\rho_w d^2}{18\mu_a}.
\end{align}
Here $f=1$ for $Re_d\leq0.01$ and $f=1+0.15Re_d^{0.687}$ otherwise,
with $Re_d<800$ enforced. This finite-slip form is regular at zero slip.
Only the two measured in-plane components are available.
The assumed still-air monodisperse references start at the lowest finite
paired mean sample on the fixed-$x$ profile, and stop at the first turning
point or upper domain. A gravity-only ballistic curve omits aerodynamic forces.
No surface launch condition, evaporation, fitted diameter or droplet tracking
is implied. A fixed-$x$ ensemble profile may decrease through lateral plume
displacement and is not a single-particle trajectory. A height-wise tracked
ensemble maximum would be a different statistic, also not a particle trajectory.

At 20 degrees C, the adopted approximate dry-air reference is
$\rho_a=1.20479\,\mathrm{kg\,m^{-3}}$; water density is evaluated from NIST
GLP10 and air viscosity from NASA's Sutherland constants (see parameter file).
Pressure 101325 Pa and equal water/air temperature are assumed. Humidity is unavailable; dry-air properties are an approximation and evaporation is neglected.
Spherical rigid-interface drag, constant diameter and one-way coupling are
model approximations; internal circulation, deformation and evaporation are omitted.

No validated diameter distribution or simultaneous gas field was found.
Diameters 5, 10, 20, 50 and 100 micrometres are sensitivity inputs.
The measured mean-field maximum is $1.686548142\,\mathrm{m\,s^{-1}}$ at
$(x,y)=(14.54665,1.284)\,\mathrm{mm}$. Quiescent-air rise distances depend
strongly on the assumed diameter; they do not validate a size or coflow.
The active Stokes-number calculation is an explicitly hypothetical jet-interaction comparison. At finite $Re_d$, $\tau_p/f$ exported along trajectories is only a
drag/slip coefficient timescale, not a linearized relaxation time.
The separate ZNMF experiment gives a hypothetical local jet-interaction comparison
$St_D(s)=\tau_p U_c(s)/D_0$, $D_0=2\,\mathrm{mm}$, not a nozzle-exit or
measured PZT Stokes number. No spatial registration is assumed.
The atomiser frequency is not assigned to the gas-flow timescale.
There is no unique measured $Re_d$ or $St$ without the missing inputs, and
a separate gas profile exceeding droplet velocities does not prove entrainment.

% Figure captions (insert each in the matching figure environment).
'''
    tex+='\nThe image-reference initial coordinate is '+f"${summary['initial_y_mm']:.6g}\\,\\mathrm{{mm}}$ with $v_0={summary['initial_mean_v_m_s']:.8g}\\,\\mathrm{{m\\,s^{{-1}}}}$.\n"
    tex+=r'All PZT comparison heights are $h=y-y_0$ with $h=0$ at this sample, not at the atomiser surface. The hypothetical jet coordinate $x_{\mathrm{ZNMF}}$ retains its separate experimental origin.'+'\n'
    tex+=f"Adopted properties are $\\rho_w={p['rho_water_kg_m3']:.7g}\\,\\mathrm{{kg\\,m^{{-3}}}}$ and $\\mu_a={p['mu_air_Pa_s']:.7g}\\,\\mathrm{{Pa\\,s}}$.\\cite{{ref1,ref2,ref3}} The drag validity follows Ref.~\\cite{{ref4}}.\n"
    tex+='\\begin{center}\n\\begin{tabular}{rrr}\nAssumed $d$ ($\\mu$m) & $\\tau_p$ (ms) & Still-air rise (mm)\\\\\n\\hline\n'
    for case in summary['cases']:
        tex+=f"{case['diameter_um']} & {case['tau_stokes_s']*1e3:.6g} & {case['end_y_mm']-summary['initial_y_mm']:.6g}"+r'\\'+'\n'
    tex+='\\end{tabular}\n\\end{center}\n% Requires graphicx; figures are relative to this section file directory.\n'
    for stem in ['velocity_theory_comparison','reynolds_vs_height','particle_relaxation_time','diameter_sensitivity','znmf_local_nozzle_stokes']:
        paragraph=text.split(stem+': ')[1].split('\n\n')[0]
        tex+='% '+stem+'\n'+r'\begin{figure}[htbp]\centering'+'\n'+r'\includegraphics[width=0.9\linewidth]{'+stem+'.pdf}\n'+r'\caption{'+paragraph.replace('_',r'\_')+'}\n'+r'\end{figure}'+'\n\n'
    tex=tex.replace('1.686548142',f"{peak['mean_v_m_s']:.10g}").replace('(14.54665,1.284)',f"({peak['x_mm']:.8g},{peak['y_mm']:.8g})")
    tex+='\\begin{thebibliography}{9}\n'
    for i,(label,url) in enumerate(SOURCES.items(),1):tex+=f'\\bibitem{{ref{i}}} {label.replace("_"," ")}. \\texttt{{\\detokenize{{{url}}}}}.\n'
    tex+='\\end{thebibliography}\n'
    (out/'methods_results.tex').write_text(tex,encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path)
    print(json.dumps(run(parser.parse_args().config),indent=2))
