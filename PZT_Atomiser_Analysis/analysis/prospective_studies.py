"""Independent dimensionless diagnostics and proposed side-jet placement study."""
import copy
import json
from pathlib import Path
import numpy as np
from .combined_flow_prediction import ROOT,PROJECT,sources_and_arrays,write_csv,write_json,sha
from .vertical_profile import extract_vertical_profile
from .droplet_dynamics import relaxation_time,reynolds
from .prescribed_gas import MeasuredGasField,integrate_droplet,trajectory_quantities
from .plotting import plt,contour_map


def save(fig,out,stem):
    for ext in ('pdf','png'):fig.savefig(out/(stem+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)


def lines(out,stem,xlabel,ylabel,curves,log=False):
    if stem.endswith('_displacement'):
        fig,axes=plt.subplots(2,1,figsize=(6.6,5.4),sharex=True,layout='constrained')
        for component,ax in enumerate(axes):
            for x,y,label,style in curves[component::2]:ax.plot(x,y,label=label.split(',')[0],lw=1.2)
            ax.set_ylabel(r'$\Delta x\ (\mathrm{mm})$' if component==0 else r'$\Delta y\ (\mathrm{mm})$')
            ax.legend(frameon=False,fontsize=8)
        axes[1].set_xlabel(xlabel);save(fig,out,stem);return
    fig,ax=plt.subplots(figsize=(6.6,4.4),layout='constrained')
    for x,y,label,style in curves:ax.plot(x,y,style,label=label,lw=1.2)
    ax.set(xlabel=xlabel,ylabel=ylabel)
    if log:ax.set_yscale('log')
    if len(curves)>1 or curves[0][2]:ax.legend(frameon=False,fontsize=8)
    save(fig,out,stem)


def dimensionless(cfg,out,pzt,D0,props):
    out=out/'Dimensionless_Numbers';out.mkdir(exist_ok=True)
    d=np.asarray(cfg['droplet_diameters_um'],float)*1e-6
    if len(d)!=1:raise ValueError('This nominal report expects one configured diameter')
    d=float(d[0]);tau=float(relaxation_time(d,props));nu=props['mu_air_Pa_s']/props['rho_air_kg_m3']
    gcfg_path=(ROOT/cfg['znmf_configuration']).resolve();gcfg=json.loads(gcfg_path.read_text());f=float(gcfg['forcing_Hz'])
    center_path=gcfg_path.parent.parent/'tables/FOV_1_2_3_centreline_and_width.csv'
    z=np.genfromtxt(center_path,names=True,delimiter=',',dtype=None,encoding='utf-8')
    good=(z['status']=='measured')&np.isfinite(z['Uc_m_per_s'])&(z['Uc_m_per_s']>0)
    uc=np.where(good,z['Uc_m_per_s'],np.nan);x=z['x_over_D0']*D0*1e3
    table=dict(x_over_D0=z['x_over_D0'],x_ZNMF_mm=x,Uc_measured_m_s=z['Uc_m_per_s'],definition_supported=good,
               Re_D_c=uc*D0/nu,Sr_D_c=f*D0/uc,Stk_D=tau*uc/D0,
               tau_Stokes_s=np.full(len(x),tau),Stk_f=np.full(len(x),tau*f),status=z['status'])
    write_csv(out/'znmf_local_centreline_numbers.csv',table)
    write_csv(out/'nominal_droplet_response.csv',dict(diameter_um=[d*1e6],tau_Stokes_s=[tau],tau_Stokes_ms=[tau*1e3],
              actuation_frequency_Hz=[f],Stk_f=[tau*f],gas_kinematic_viscosity_m2_s=[nu],D0_m=[D0]))
    rp=reynolds(d,np.stack([pzt['mean_u'],pzt['mean_v']],axis=-1),[0,0],props)
    write_csv(out/'pzt_zero_background_reference_field.csv',dict(x_mm=pzt['x_mm'].ravel(),y_mm=pzt['y_mm'].ravel(),
        mean_u_m_s=pzt['mean_u'].ravel(),mean_v_m_s=pzt['mean_v'].ravel(),N=pzt['N'].ravel(),Re_p_zero_gas_reference=rp.ravel()))
    profile,peak=extract_vertical_profile(ROOT/cfg['pzt_statistics']);i,j=profile['native_row'],profile['native_column']
    write_csv(out/'pzt_fixed_line_zero_background_reference.csv',dict(**profile,mean_u=pzt['mean_u'][i,j],Re_p_zero_gas_reference=rp[i,j]))
    contour_map(pzt['x_mm'],pzt['y_mm'],rp,label=r'$Re_p^{\mathrm{ref}}\ (-)$',limits=(0,float(np.nanmax(rp))),
        output=out/'pzt_Re_p_zero_background_reference',xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$',
        data_output=out/'pzt_Re_p_zero_background_reference.npz')
    lines(out,'pzt_profile_Re_p_zero_background_reference',r'$y\ (\mathrm{mm})$',r'$Re_p^{\mathrm{ref}}\ (-)$',
          [(profile['y_mm'],rp[i,j],rf'Fixed $x={peak["x_mm"]:.5f}$ mm','r-')])
    for key,label in [('Re_D_c',r'$Re_{D,c}(x)\ (-)$'),('Sr_D_c',r'$Sr_{D,c}(x)\ (-)$'),('Stk_D',r'$Stk_D(x)\ (-)$')]:
        lines(out,'znmf_'+key,r'$x_{\mathrm{ZNMF}}\ (\mathrm{mm})$',label,[(x,table[key],'Local measured centreline','-')])
    lines(out,'droplet_Stk_f',r'$x_{\mathrm{ZNMF}}\ (\mathrm{mm})$',r'$Stk_f=\tau_p f\ (-)$',
          [(x,np.where(good,tau*f,np.nan),'Frequency response diagnostic','-')])
    summary=dict(diameter_um=d*1e6,tau_Stokes_s=tau,actuation_frequency_Hz=f,acquisition_frequency_Hz_not_used=gcfg['acquisition_Hz'],
                 Stk_f=tau*f,nu_air_m2_s=nu,D0_m=D0,centreline_valid_rows=int(good.sum()),centreline_total_rows=len(x),
                 source=str(center_path),source_sha256=sha(center_path),
                 formulas=dict(Re_p='rho_air*d*norm(measured mean droplet vector)/mu_air under explicit zero-background reference',
                               Re_D_c='Uc*D0/nu',Sr_D_c='f*D0/Uc',Stk_f='tau_Stokes*f, no 2*pi',Stk_D='tau_Stokes*Uc/D0'))
    write_json(out/'parameters_and_validation.json',summary)
    np.testing.assert_allclose(table['Re_D_c'][good]*table['Sr_D_c'][good],f*D0**2/nu,rtol=1e-14)
    np.testing.assert_allclose(table['Stk_D'][good]*table['Sr_D_c'][good],tau*f,rtol=1e-14)
    assert np.isnan(table['Re_D_c'][~good]).all()
    tex=rf'''\section{{Dimensionless diagnostics from separate experiments}}
The configured $d=10\,\mu\mathrm{{m}}$ is a nominal model estimate, motivated by
$1.89$ times a nominal $5\,\mu\mathrm{{m}}$ aperture under an equal-jet-diameter
assumption, not measured sizing. The Stokes-limit response time is
$\tau_p=\rho_w d^2/(18\mu_a)={tau*1e3:.6g}\,\mathrm{{ms}}$.
PZT-only $Re_p=\rho_a d\lVert\overline{{\boldsymbol w}}_d\rVert/\mu_a$ uses an
explicit zero-background-gas reference. It is conditional, not a measurement of
ambient gas or a claim that the PZT experiment had stationary surrounding air.
Signed measured components are retained; magnitude is used only in defining Re.
This Re of a mean vector is not an ensemble mean of instantaneous Re.

For the separately measured ZNMF jet, $Re_{{D,c}}(x)=U_c(x)D_0/\nu_a$ and
$Sr_{{D,c}}(x)=fD_0/U_c(x)$ are local centreline-based quantities with
$D_0=2\,$mm and documented actuation $f={f:g}\,$Hz, not the
{gcfg['acquisition_Hz']:g} Hz acquisition rate or PIV pulse frequency.
Only finite, positive, measured $U_c$ is used; all other rows remain gaps.
No nozzle-exit velocity, $Re_0$ or $Sr_0$ is inferred.

Droplet response diagnostics are $Stk_f=\tau_p f={tau*f:.7g}$ (without $2\pi$)
and $Stk_D(x)=\tau_p U_c(x)/D_0$. Stokes number Stk is distinct from
Strouhal number Sr. The response time is a Stokes-limit baseline, not an exact
finite-Re relaxation time. These separate-experiment diagnostics are not
integrated trajectories, evidence of capture, or measurements of the combined flow.
The design study separately evaluates evolving vector slip for trajectory Re.

\paragraph{{Figure captions.}}
PZT maps/profiles: conditional $Re_p$ from signed measured mean components and
the explicit zero-background-gas reference, at the nominal assumed diameter.
ZNMF local-number figures: measured-centreline definitions retaining unsupported
gaps, with frequency from the actuation configuration. Frequency-response figure:
constant $Stk_f$, shown only on supported measured-centreline rows.
'''
    (out/'methods_results.tex').write_text(tex,encoding='utf-8')
    stems=['pzt_Re_p_zero_background_reference','pzt_profile_Re_p_zero_background_reference','znmf_Re_D_c','znmf_Sr_D_c','znmf_Stk_D','droplet_Stk_f']
    write_json(out/'manifest.json',dict(figures=stems,CSVs=[p.name for p in out.glob('*.csv')],geometry_independent=True))
    return summary


def design(cfg,out,pzt,jet,xj,yj,D0,props):
    out=out/'Design_Study';out.mkdir(exist_ok=True)
    settings=cfg['design_study'];controls=copy.deepcopy(cfg);controls['solver']['maximum_time_s']=settings['maximum_time_s']
    d=float(cfg['droplet_diameters_um'][0])*1e-6
    profile,peak=extract_vertical_profile(ROOT/cfg['pzt_statistics']);anchor_x=peak['x_mm'];col=peak['native_column']
    with np.load(ROOT/cfg['znmf_field']) as z:
        valid=z['centre_valid']&np.isfinite(z['Uc'])&(z['Uc']>0)&np.isfinite(z['yc_measured'])
        supported_columns=np.flatnonzero(valid);centres=z['yc_measured'].copy()
    ylo,yhi=float(np.nanmin(pzt['y_mm'])),float(np.nanmax(pzt['y_mm']))
    gas_origin=float(json.loads((ROOT/cfg['znmf_nozzle_reference']).read_text())['geometric_nozzle_centre_registered_D0'])
    cases=[];placements=[];convergence=[];manifest=[]
    from matplotlib.colors import Normalize
    for ih,hfraction in enumerate(settings['plume_height_fractions']):
        target=ylo+hfraction*(yhi-ylo);row=int(np.argmin(np.abs(pzt['y_mm'][:,col]-target)));height=float(pzt['y_mm'][row,col])
        for ix,xfraction in enumerate(settings['jet_downstream_support_fractions']):
            index=int(supported_columns[round(xfraction*(len(supported_columns)-1))]);distance=float(xj[index]*1e3)
            local_centre_y=float((centres[index]-gas_origin)*D0*1e3)
            nozzle=np.array([anchor_x-distance,height-local_centre_y])
            field=MeasuredGasField(xj,yj,jet['U'],jet['V'],nozzle*1e-3,0)
            placement=f'H{ih+1}_L{ix+1}';placement_cases=[]
            for offset in settings['release_column_offsets']:
                j=col+offset;initial=np.array([pzt['x_mm'][row,j]*1e-3,pzt['y_mm'][row,j]*1e-3,pzt['mean_u'][row,j],pzt['mean_v'][row,j]])
                if not np.isfinite(initial).all():raise ValueError('Selected design release lacks finite measured PZT support')
                field.sample(initial[:2]);ident=f'{placement}_C{j}'
                tr,reason=integrate_droplet(field,initial,d,props,controls)
                fine,fine_reason=integrate_droplet(field,initial,d,props,controls,.5)
                end=min(tr[-1,0],fine[-1,0]);times=np.linspace(0,end,201)
                errors=[float(np.max(np.abs(np.interp(times,tr[:,0],tr[:,k])-np.interp(times,fine[:,0],fine[:,k])))) for k in range(1,5)]
                if max(errors[:2])>cfg['solver']['position_convergence_tolerance_m'] or max(errors[2:])>cfg['solver']['velocity_convergence_tolerance_m_s']:raise AssertionError((ident,'step convergence',errors))
                if len(tr)<2 or tr[-1,0]<=0:raise AssertionError('No successful integration')
                q=trajectory_quantities(field,tr,d,props,D0,cfg)
                q['horizontal_deflection_mm']=q['x_mm']-q['x_mm'][0];q['vertical_displacement_mm']=q['y_mm']-q['y_mm'][0]
                write_csv(out/(ident+'.csv'),q)
                record=dict(case=ident,placement=placement,proposed_geometry=True,plume_centreline_height_mm=height,
                    nozzle_x_mm=float(nozzle[0]),nozzle_y_mm=float(nozzle[1]),nozzle_to_plume_distance_mm=distance,
                    release_row=row,release_column=j,initial_x_mm=float(initial[0]*1e3),initial_y_mm=float(initial[1]*1e3),
                    initial_u_m_s=float(initial[2]),initial_v_m_s=float(initial[3]),diameter_um=d*1e6,
                    termination=reason,end_time_s=float(tr[-1,0]),horizontal_deflection_mm=float(q['horizontal_deflection_mm'][-1]),
                    vertical_displacement_mm=float(q['vertical_displacement_mm'][-1]),max_particle_Re=float(q['Re_in_plane'].max()))
                cases.append(record);placement_cases.append((ident,q));convergence.append(dict(case=ident,half_step_termination=fine_reason,max_errors_x_y_u_v=errors))
                print('Integrated',ident,reason,flush=True)
            xx,yy,ug,vg=field.common_mesh();extent=float(np.nanmax(np.abs(ug)))
            fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
            m=ax.contourf(xx*1e3,yy*1e3,np.ma.masked_invalid(ug),levels=np.linspace(-extent,extent,41),norm=Normalize(-extent,extent),cmap='RdBu_r',corner_mask=False)
            fig.colorbar(m,ax=ax,label=r'$\overline{u}_g\ (\mathrm{m\,s^{-1}})$',shrink=.75)
            for ident,q in placement_cases:ax.plot(q['x_mm'],q['y_mm'],lw=1.5,label=ident)
            ax.set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$');ax.set_aspect('equal');ax.legend(frameon=False,fontsize=8)
            stem=placement+'_trajectories';save(fig,out,stem);manifest.append(stem)
            for stem,cols,label in [('velocity',['droplet_u_m_s','droplet_v_m_s'],r'$u_d,\ v_d\ (\mathrm{m\,s^{-1}})$'),
                                   ('displacement',['horizontal_deflection_mm','vertical_displacement_mm'],r'$\Delta x,\ \Delta y\ (\mathrm{mm})$'),
                                   ('particle_Re',['Re_in_plane'],r'$Re_p\ (-)$'),('supporting_Stokes',['St_D_local'],r'$Stk_{D,\mathrm{local}}\ (-)$')]:
                curves=[]
                for ident,q in placement_cases:
                    for k,key in enumerate(cols):curves.append((q['t_s']*1e3,q[key],ident+(', '+('horizontal' if k==0 else 'vertical') if len(cols)>1 else ''),'-' if k==0 else '--'))
                name=placement+'_'+stem;lines(out,name,r'$t\ (\mathrm{ms})$',label,curves);manifest.append(name)
            relevant=[c for c in cases if c['placement']==placement]
            placements.append(dict(placement=placement,proposed_geometry=True,centreline_height_mm=height,nozzle_to_plume_distance_mm=distance,
                nozzle_x_mm=float(nozzle[0]),nozzle_y_mm=float(nozzle[1]),nominal_nozzle_axis_height_mm=float(nozzle[1]),
                integrated_cases=len(relevant),termination_reasons='; '.join(c['termination'] for c in relevant),
                min_horizontal_deflection_mm=min(c['horizontal_deflection_mm'] for c in relevant),max_horizontal_deflection_mm=max(c['horizontal_deflection_mm'] for c in relevant),
                min_vertical_displacement_mm=min(c['vertical_displacement_mm'] for c in relevant),max_vertical_displacement_mm=max(c['vertical_displacement_mm'] for c in relevant)))
    for name,records in [('case_comparison.csv',cases),('placement_comparison.csv',placements)]:write_csv(out/name,{k:[r[k] for r in records] for k in records[0]})
    write_json(out/'convergence.json',convergence)
    write_json(out/'manifest.json',dict(successfully_integrated_cases=len(cases),placement_count=len(placements),figures=manifest,
        trajectory_CSVs=[c['case']+'.csv' for c in cases],settings=settings,measured_PZT_y_range_mm=[ylo,yhi],
        measured_ZNMF_downstream_range_mm=[float(xj[supported_columns[0]]*1e3),float(xj[supported_columns[-1]]*1e3)],
        selection='Nearest native PZT row at configured fractions of full measured height; gas downstream stations at fractions of positive measured centreline station indices. Align that measured gas centreline to PZT anchor, not the nominal nozzle axis. Three adjacent native PZT columns initialise both signed mean components.'))
    (out/'methods_results.tex').write_text(r'''\section{Prospective horizontal side-jet design study}
Relative placement is a design variable, not measured geometry. Nine proposed
arrangements combine three heights within the measured PZT domain with three
nozzle-to-plume distances within the measured positive-centreline ZNMF support.
Selected measured gas-centreline positions are aligned to the existing PZT
fixed-line anchor; nozzle-axis heights differ by the measured gas-centreline offset.
Actual native indices, coordinates, nozzle placements and initial signed velocity
components are exported in the comparison tables. Three adjacent supported native
PZT columns per placement supply 27 representative releases at nominal $d=10\,\mu$m.
These are conditional initialisations in the interaction region, not simulated
unmeasured travel from the atomiser surface.

The separately measured two-component ZNMF mean gas field is prescribed after
rotation/translation of positions and vectors. PZT mean components initialise
particles once and are never added to the gas or repeatedly imposed thereafter.
Finite-Re spherical drag acts on evolving vector slip, with gravity and buoyancy;
constant diameter and no evaporation are assumed. This one-way mean-field model
neglects plume-induced gas changes and jet unsteadiness. It is not experimentally
validated for combined flow. Source cells require four finite vectors; no gas
extrapolation is used and trajectories stop at unsupported boundaries or the
configured 30 ms time horizon. Termination reasons are tabulated. Half-step
comparisons verify position and velocity convergence.

Trajectory Reynolds number uses instantaneous model gas-droplet vector slip.
Supporting $Stk_{D,local}=\tau_{p,Stokes}\lVert\boldsymbol w_g\rVert/D_0$
uses the local nozzle-length gas timescale, not nozzle-exit conditions or the
atomiser frequency. It is not a capture criterion. Deflections and displacements
refer to each release point. No collection target is defined, so no capture
efficiency or optimal arrangement is claimed.

\paragraph{Figure captions.} Trajectory panels: proposed placements over
separately measured mean horizontal gas-velocity contours; curves are integrated
conditional predictions. Velocity panels: horizontal (solid) and vertical (dashed) components for the same releases. Displacement panels use separate horizontal and vertical ordinates. Particle-Re panels use
evolving vector slip; supporting Stokes panels use the explicitly local gas timescale.
''',encoding='utf-8')
    return dict(successfully_integrated_cases=len(cases),placement_count=len(placements))


def run_dimensionless():
    """Generate dimensionless outputs with no placement or trajectory prerequisites."""
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    if cfg.get('co_directed',{}).get('enabled'):
        from .co_directed import context
        from .co_directed_report import dimensionless as upward_numbers
        c,p,j,D,s,props,active,profile,peak,field,alignment=context()
        return upward_numbers(active,c,props,D,alignment)[0]
    pzt,jet,xj,yj,D0,source=sources_and_arrays(cfg)
    props_path=out/'fluid_properties.json';props=json.loads(props_path.read_text())
    result=dimensionless(cfg,out,pzt,D0,props)
    write_json(out/'Dimensionless_Numbers/configuration_used.json',cfg)
    write_json(out/'Dimensionless_Numbers/source_manifest.json',dict(**source,saved_fluid_properties_path=str(props_path),saved_fluid_properties_sha256=sha(props_path)))
    return result


def run():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    if cfg.get('co_directed',{}).get('enabled'):
        from .co_directed import run as upward_run
        return upward_run()
    pzt,jet,xj,yj,D0,source=sources_and_arrays(cfg)
    props_path=out/'fluid_properties.json';props=json.loads(props_path.read_text())
    before={s['path']:s['sha256'] for s in source['sources']}
    numbers=dimensionless(cfg,out,pzt,D0,props)  # No placement gate.
    if not cfg.get('design_study',{}).get('enabled'):return dict(design=None,dimensionless=numbers)
    study=design(cfg,out,pzt,jet,xj,yj,D0,props)
    for path,expected in before.items():assert sha(path)==expected,path
    write_json(out/'prospective_studies_status.json',dict(design=study,dimensionless=numbers,source_files_unchanged=True))
    status=json.loads((out/'prediction_status.json').read_text())
    status['fixed_user_placement_case_count']=status.get('prediction_case_count',0)
    status['prediction_case_count']=study['successfully_integrated_cases']
    status['proposed_design_placement_count']=study['placement_count']
    status['status']='prospective_design_study_complete'
    status['pending_fixed_geometry_inputs']=status.pop('missing_inputs',[])
    status['missing_inputs']=[]
    status['geometry_status']='Proposed design variables selected from measured support; not measured combined geometry'
    write_json(out/'prediction_status.json',status)
    report=json.loads((out/'report_manifest.json').read_text())
    report['prospective_design_manifest']='Design_Study/manifest.json';report['dimensionless_manifest']='Dimensionless_Numbers/manifest.json'
    write_json(out/'report_manifest.json',report)
    for folder in ['Dimensionless_Numbers','Design_Study']:
        write_json(out/folder/'configuration_used.json',cfg);write_json(out/folder/'source_manifest.json',dict(**source,saved_fluid_properties_path=str(props_path),saved_fluid_properties_sha256=sha(props_path)))
    refresh_reports()
    return dict(design=study,dimensionless=numbers)


def refresh_reports():
    """Summarise completed CSV results without re-integrating trajectories."""
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    status=json.loads((out/'prediction_status.json').read_text());sources=json.loads((out/'source_manifest.json').read_text())
    from .combined_flow_report import write_report
    write_report(out,cfg,sources,status,json.loads((out/'report_manifest.json').read_text())['figures'])
    placement=np.genfromtxt(out/'Design_Study/placement_comparison.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    heights=np.unique(placement['centreline_height_mm']);distances=np.unique(placement['nozzle_to_plume_distance_mm'])
    numbers=json.loads((out/'Dimensionless_Numbers/parameters_and_validation.json').read_text())
    study=out/'Design_Study'
    (study/'README.md').write_text(f'''# Prospective horizontal side-jet placements

Successfully integrated: {status['prediction_case_count']} representative droplets,
across {len(placement)} PROPOSED arrangements. This is not measured combined geometry.
Jet centreline heights: {heights.tolist()} mm, selected at nearest native PZT rows
to 25%, 50% and 75% of the measured height range. Nozzle-to-plume distances:
{distances.tolist()} mm, selected at 20%, 50% and 80% of positive measured gas
centreline station indices. The measured local gas centreline, not its nominal
nozzle axis, is aligned at the PZT fixed-line anchor. Three adjacent native columns
41, 42, 43 initialise both measured signed components at each selected height.
No travel from an unmeasured atomiser surface is simulated.

All 27 cases reach the downstream measured-gas boundary before the 30 ms limit.
Endpoints are support limits, not droplet stopping or capture. Their displacements
occur over different elapsed times and must not be ranked as collection performance.
No collection target is specified and no optimum or capture efficiency is claimed.

Files: placement_comparison.csv (nine arrangements), case_comparison.csv (27 cases,
initial vectors, displacement and termination), H*_L*_C*.csv (27 full histories).
For each H1..H3 / L1..L3 pair, PDF and 600-dpi PNG figures are:
*_trajectories, *_velocity, *_displacement, *_particle_Re, *_supporting_Stokes.
Displacement components have separate ordinates; velocities retain signed values.
CSV histories also contain both gas components, vector slip, Re and defined response
timescales. Supporting Stokes uses local D0 / gas-vector magnitude, not centreline
Uc or a nozzle-exit timescale. See methods_results.tex for the prescribed-field
model approximation, source_manifest.json for sources, manifest.json for exact
selection/configuration and convergence.json for all half-step comparisons.

Regenerate with the single command in the parent README. To refresh notebook
display from existing computed CSV/figures without new integration, use
update_combined_flow_notebook --reuse-saved-studies. This is a display refresh only.
''',encoding='utf-8')
    dim=out/'Dimensionless_Numbers'
    (dim/'README.md').write_text(f'''# Geometry-independent dimensionless diagnostics

Nominal model diameter 10 um; tau_Stokes={numbers['tau_Stokes_s']:.12g} s;
documented ZNMF actuation f=120 Hz; Stk_f={numbers['Stk_f']:.12g} (no 2*pi).
Acquisition 10 Hz and atomiser 190 kHz are not used as the jet actuation frequency.

CSV files:
- nominal_droplet_response.csv: diameter, Stokes-limit time, nu_air, D0, f, Stk_f.
- pzt_zero_background_reference_field.csv: signed measured mean u/v, coordinates,
  N and conditional particle Re under an explicit zero-background-gas reference.
- pzt_fixed_line_zero_background_reference.csv: original fixed-line profile,
  native indices and the same reference Re. Ambient PZT air was not measured.
- znmf_local_centreline_numbers.csv: all 192 original rows, measured Uc, support
  flags, local Re_D,c, Sr_D,c, Stk_D and Stk_f. Ten unsupported rows remain gaps.

PDF/600-dpi PNG stems:
pzt_Re_p_zero_background_reference; pzt_profile_Re_p_zero_background_reference;
znmf_Re_D_c; znmf_Sr_D_c; znmf_Stk_D; droplet_Stk_f.

Re_D,c=Uc*D0/nu; Sr_D,c=f*D0/Uc; Stk_D=tau_Stokes*Uc/D0. These are local
centreline-based diagnostics, never Re0/Sr0 or measured combined-flow trajectories.
Stokes Stk is distinct from Strouhal Sr. See the single methods_results.tex section,
parameters_and_validation.json and source_manifest.json for values and provenance.

Run independently of ALL placement settings:
.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.prospective_studies --dimensionless-only
''',encoding='utf-8')
    methods=study/'methods_results.tex'
    text=methods.read_text(encoding='utf-8').split('% Actual completed study results')[0]
    text+=r'''
% Actual completed study results
\paragraph{Completed results.} All 27 integrations succeeded across the nine
proposed placements. All terminate at the downstream measured-gas boundary,
not at a physical stopping condition. Endpoint deflections compare different
supported elapsed times and are not a ranking of collection performance.
The drag multiplier is unity for $Re_p\leq0.01$ and
$1+0.15Re_p^{0.687}$ otherwise, with $Re_p<800$ enforced.
'''
    methods.write_text(text,encoding='utf-8')


def notebook_cells():
    import nbformat as nb
    text='''## Prospective placement study and independent dimensionless diagnostics

Placement is a design variable: nine proposed arrangements, with three native
measured PZT initial vectors each. No measured combined geometry or capture
efficiency is implied. Independent dimensionless outputs use nominal 10 um,
saved fluid properties, full authoritative centreline support and documented
120 Hz jet actuation. Stokes Stk and Strouhal Sr are distinct.
'''
    code='''study_root=output/'Combined_Flow_Prediction'
print((study_root/'prospective_studies_status.json').read_text())
for name in ['pzt_profile_Re_p_zero_background_reference','znmf_Re_D_c','znmf_Sr_D_c','znmf_Stk_D','droplet_Stk_f']:
    display(Image(filename=str(study_root/'Dimensionless_Numbers'/(name+'.png')),width=650))
design_manifest=json.loads((study_root/'Design_Study/manifest.json').read_text())
trajectory_stems=[s for s in design_manifest['figures'] if s.endswith('_trajectories')]
for stem in trajectory_stems:
    display(Image(filename=str(study_root/'Design_Study'/(stem+'.png')),width=650))
middle=trajectory_stems[len(trajectory_stems)//2].removesuffix('_trajectories')
for suffix in ['velocity','displacement','particle_Re','supporting_Stokes']:
    display(Image(filename=str(study_root/'Design_Study'/(middle+'_'+suffix+'.png')),width=650))
'''
    cells=[nb.v4.new_markdown_cell(text),nb.v4.new_code_cell(code)]
    for cell in cells:cell.metadata['prospective_studies']=True
    return cells


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--dimensionless-only',action='store_true')
    args=parser.parse_args();print(json.dumps(run_dimensionless() if args.dimensionless_only else run(),indent=2))
