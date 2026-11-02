"""PZT-only baseline and input-gated prediction in a separate measured side jet."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from .plotting import plt,contour_map,save_figure
from .vertical_profile import extract_vertical_profile
from .droplet_dynamics import fluid_properties,reynolds,relaxation_time
from .prescribed_gas import MeasuredGasField,integrate_droplet,trajectory_quantities,UnsupportedPosition

ROOT=Path(__file__).resolve().parents[2]
PROJECT=ROOT/'PZT_Atomiser_Analysis'


def sha(path):
    path=str(Path(path).resolve())
    if os.name=='nt' and not path.startswith('\\\\?\\'):path='\\\\?\\'+path
    with open(path,'rb') as f:
        h=hashlib.sha256()
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def write_csv(path,columns):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(columns);w.writerows(zip(*columns.values()))


def sources_and_arrays(cfg):
    paths={k:(ROOT/cfg[k]).resolve() for k in ('pzt_statistics','pzt_configuration','znmf_field','znmf_configuration','znmf_nozzle_reference')}
    pcfg=json.loads(paths['pzt_configuration'].read_text());gcfg=json.loads(paths['znmf_configuration'].read_text())
    ref=json.loads(paths['znmf_nozzle_reference'].read_text())
    np.testing.assert_allclose([pcfg['coordinate_scale'],pcfg['velocity_scale'],pcfg['confirmed_pulse_separation_seconds']],[.00535,.214,25e-6],rtol=0,atol=1e-15)
    with np.load(paths['pzt_statistics']) as z:pzt={k:z[k].copy() for k in z.files}
    with np.load(paths['znmf_field']) as z:jet={k:z[k].copy() for k in ('x','y','U','V')}
    D0=gcfg['diameter_mm']*1e-3;nozzle_y=ref['geometric_nozzle_centre_registered_D0']
    local_x=jet['x']*D0;local_y=(jet['y']-nozzle_y)*D0
    if not ref['geometric_nozzle_centre_established']:raise ValueError('Jet nozzle-reference evidence unavailable')
    manifest=dict(sources=[dict(role=k,path=str(v),sha256=sha(v)) for k,v in paths.items()],
        experiments=dict(PZT='Atomiser only, 190 kHz/9 Vpp; no ZNMF jet; measured droplet-phase velocities; ambient gas unmeasured',
                         ZNMF='Air jet only without PZT droplets; mean two-component gas velocities',combined='Not conducted; predictions are unvalidated'),
        pzt=dict(coordinates='Saved x_mm,y_mm retained in image-reference frame, x right/y up',velocities='Saved mean_u,mean_v and frame u,v read directly in m/s',
                 scaling='No velocity recalibration; coordinates converted mm to m only for integration',population=pcfg['working_assumptions'][-1],
                 finite_paired_vectors=int(np.sum(np.isfinite(pzt['mean_u'])&np.isfinite(pzt['mean_v'])))),
        znmf=dict(coordinates='Local Xj=x*D0; Yj=(y-6.868)*D0, from registered dimensionless arrays; no use of y_plot',
                  nozzle_y_registered_D0=float(nozzle_y),D0_m=D0,velocity='U,V already m/s; downstream sign reversal already applied upstream, not applied again',
                  masks='Both stored U,V finite; preserve upstream reconstruction support; bilinear evaluation only in cells with four finite paired vectors',
                  finite_paired_vectors=int(np.sum(np.isfinite(jet['U'])&np.isfinite(jet['V']))),
                  local_x_range_mm=(local_x[[0,-1]]*1e3).tolist(),local_y_range_mm=(local_y[[0,-1]]*1e3).tolist(),
                  reconstruction='Authoritative three-FOV reconstructed mean, including its existing overlap registration/blending; not reprocessed'),
        placement=dict(common_frame=cfg['common_frame'],transformation='q_common=R(theta)*q_jet_local+nozzle_position; w_common=R(theta)*(U,V)',
                       pzt_transform='Identity: common frame retains final PZT right/up axes and image-reference origin',
                       geometry=cfg.get('geometry',{}),warning='Origins of the two experiments are not equated'))
    return pzt,jet,local_x,local_y,D0,manifest


def missing_inputs(cfg):
    missing=[]
    if not cfg['droplet_diameters_um']:missing.append('Droplet diameter(s) in micrometres, with sizing evidence or an explicit parameter-study assumption')
    elif not cfg.get('diameter_source_or_assumption'):missing.append('State whether the configured diameter values are measured or assumed, with their source')
    geometry=cfg['geometry']
    if geometry['jet_nozzle_position_in_pzt_frame_mm'] is None or geometry['jet_rotation_degrees'] is None or not geometry.get('placement_source'):
        missing.append('Proposed ZNMF nozzle-centre (x,y) in the PZT image-reference frame and jet direction/rotation (including which side it enters), with placement description')
    if not cfg['release_locations']:missing.append('Droplet entry/release location(s): native PZT row/column or measured (x,y) identifying those nodes, and ensemble-mean or saved-frame initial velocity selection')
    return missing


def _record(manifest,stem,kind,caption):
    manifest.append(dict(stem=stem,kind=kind,exports=[f'figures/{stem}.{e}' for e in ('pdf','svg','png')],caption=caption))


def baseline(out,pzt,jet,xj,yj,source,manifest):
    data=out/'data';figures=out/'figures'
    data.mkdir(exist_ok=True);figures.mkdir(exist_ok=True)
    for name in ('mean_u','mean_v'):
        value=pzt[name];extent=float(np.nanmax(np.abs(value)))
        label=r'$\overline{u}\ (\mathrm{m\,s^{-1}})$' if name=='mean_u' else r'$\overline{v}\ (\mathrm{m\,s^{-1}})$'
        stem='pzt_only_'+name
        contour_map(pzt['x_mm'],pzt['y_mm'],value,label=label,limits=(-extent,extent),signed=True,
                    output=figures/stem,xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$',data_output=data/(stem+'.npz'))
        _record(manifest,stem,'PZT-only measured droplet baseline',
                f'PZT-only ensemble {name} in the preserved image-reference frame. The ZNMF jet was absent; surrounding-air velocity was not measured. All signed values and existing finite-data support are retained. No gas velocity is inferred.')
    write_csv(data/'pzt_only_mean_velocity_field.csv',dict(x_mm=pzt['x_mm'].ravel(),y_mm=pzt['y_mm'].ravel(),
        mean_u_m_s=pzt['mean_u'].ravel(),mean_v_m_s=pzt['mean_v'].ravel(),N=pzt['N'].ravel(),
        finite_pair=(np.isfinite(pzt['mean_u'])&np.isfinite(pzt['mean_v'])).ravel()))
    profile,peak=extract_vertical_profile(source)
    rows=profile['native_row'];cols=profile['native_column'];u=pzt['mean_u'][rows,cols]
    write_csv(data/'pzt_only_fixed_line_profile.csv',dict(**profile,mean_u=u))
    np.savez_compressed(data/'pzt_only_fixed_line_profile.npz',**profile,mean_u=u)
    for component,values,color in [('u',u,'tab:blue'),('v',profile['mean_v'],'#d62728')]:
        stem='pzt_only_mean_'+component+'_vertical_profile';fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
        ax.plot(profile['y_mm'],values,color=color,lw=1.3,label=rf'$x={peak["x_mm"]:.5f}\,\mathrm{{mm}}$')
        ax.set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=rf'$\overline{{{component}}}\ (\mathrm{{m\,s^{{-1}}}})$');ax.legend(frameon=False)
        save_figure(fig,figures/stem,dict(y_mm=profile['y_mm'],mean_velocity=values),data_output=data/(stem+'.npz'))
        _record(manifest,stem,'PZT-only measured droplet baseline',
                f'Fixed-x ensemble-mean droplet {component} profile through the existing global mean-v maximum, x={peak["x_mm"]:.5f} mm. One native column is retained over all measured y values, preserving signs and gaps. This is not an individual droplet trajectory or tracked plume centreline. The surrounding-air velocity is unknown.')
    # Native gas-only views need no proposed alignment and are never overlaid on PZT.
    xx,yy=np.meshgrid(xj*1e3,yj*1e3)
    for key in ('U','V'):
        values=jet[key];extent=float(np.nanmax(np.abs(values)));stem='znmf_only_mean_'+key
        contour_map(xx,yy,values,label=rf'$\overline{{{key}}}_g\ (\mathrm{{m\,s^{{-1}}}})$',limits=(-extent,extent),signed=True,
                    output=figures/stem,xlabel=r'$X_j\ (\mathrm{mm})$',ylabel=r'$Y_j\ (\mathrm{mm})$',data_output=data/(stem+'.npz'),figsize=(7.2,3.4))
        _record(manifest,stem,'ZNMF-only measured gas input',
                'Separate ZNMF-only mean gas component in its own nozzle-centred local frame. Xj is downstream and Yj transverse/upward. Existing velocity signs and reconstruction masks are retained. No placement in the PZT frame has been assumed.')
    np.savez_compressed(data/'znmf_only_local_mean_field.npz',X_j_m=xj,Y_j_m=yj,U_m_s=jet['U'],V_m_s=jet['V'],finite_pair=np.isfinite(jet['U'])&np.isfinite(jet['V']))
    write_csv(data/'znmf_only_local_mean_field.csv',dict(X_j_mm=xx.ravel(),Y_j_mm=yy.ravel(),
        x_registered_over_D0=np.broadcast_to(jet['x'],jet['U'].shape).ravel(),
        y_registered_over_D0=np.broadcast_to(jet['y'][:,None],jet['U'].shape).ravel(),
        U_m_s=jet['U'].ravel(),V_m_s=jet['V'].ravel(),finite_pair=(np.isfinite(jet['U'])&np.isfinite(jet['V'])).ravel()))
    return peak


def release_state(release,pzt,statistics_path):
    i,j=int(release['row']),int(release['column'])
    if not (0<=i<pzt['mean_v'].shape[0] and 0<=j<pzt['mean_v'].shape[1]):raise ValueError('Release native index outside PZT field')
    kind=release.get('velocity_source','ensemble_mean')
    if kind=='ensemble_mean':u,v=pzt['mean_u'][i,j],pzt['mean_v'][i,j]
    elif kind.startswith('frame:'):
        key=kind.split(':',1)[1]
        if not key.isdigit():raise ValueError('Saved frame key must be numeric')
        with np.load(statistics_path.parent/f'frame_{key}.npz') as z:u,v=z['u'][i,j],z['v'][i,j]
    else:raise ValueError('Use ensemble_mean or frame:<saved key>')
    if not (np.isfinite(u) and np.isfinite(v)):raise ValueError('Release velocity pair is not finite')
    return np.array([pzt['x_mm'][i,j]*1e-3,pzt['y_mm'][i,j]*1e-3,u,v])


def prediction_figures(out,field,cases,manifest):
    figures=out/'figures';data=out/'data'
    xx,yy,ug,vg=field.common_mesh();extent=float(np.nanmax(np.abs(ug)))
    fig,ax=plt.subplots(figsize=(7.2,4.8),layout='constrained')
    levels=np.linspace(-extent,extent,41)
    from matplotlib.colors import Normalize
    c=ax.contourf(xx*1e3,yy*1e3,np.ma.masked_invalid(ug),levels=levels,cmap='RdBu_r',norm=Normalize(-extent,extent),corner_mask=False)
    fig.colorbar(c,ax=ax,label=r'$\overline{u}_g\ (\mathrm{m\,s^{-1}})$')
    arrays={}
    for case in cases:
        q=case['quantities'];label=case['label'];ax.plot(q['x_mm'],q['y_mm'],label=label,lw=1.3)
        arrays[case['id']+'_x_mm']=q['x_mm'];arrays[case['id']+'_y_mm']=q['y_mm']
    ax.set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$');ax.set_aspect('equal');ax.legend(fontsize=8,frameon=False)
    save_figure(fig,figures/'predicted_trajectories',dict(**arrays,gas_x_m=xx,gas_y_m=yy,gas_u_m_s=ug,gas_v_m_s=vg),data_output=data/'predicted_trajectories.npz')
    _record(manifest,'predicted_trajectories','Unvalidated combined-arrangement prediction','Representative droplet predictions over the prescribed separately measured ZNMF mean gas field in explicitly configured placement. PZT velocities initialise the particles only. The mean field is not a measured combined-flow gas field. Curves stop at unsupported measured gas boundaries/gaps, drag-law limits or configured time limit.')
    specifications=[('predicted_velocity_components',['droplet_u_m_s','droplet_v_m_s'],r'$u_d,\ v_d\ (\mathrm{m\,s^{-1}})$'),
                    ('predicted_slip_components',['slip_u_m_s','slip_v_m_s','slip_magnitude_m_s'],r'$u_{a}-u_d,\ v_a-v_d,\ |\boldsymbol w_a-\boldsymbol w_d|\ (\mathrm{m\,s^{-1}})$'),
                    ('predicted_particle_reynolds',['Re_in_plane'],r'$Re_d\ (-)$'),
                    ('predicted_local_stokes',['St_D_local'],r'$St_{D,\mathrm{local}}\ (-)$')]
    for stem,columns,ylabel in specifications:
        fig,ax=plt.subplots(figsize=(7.2,4.8),layout='constrained');arrays={}
        for case in cases:
            q=case['quantities']
            for column,style in zip(columns,('-','--',':')):
                ax.plot(q['t_s']*1e3,q[column],ls=style,label=case['label']+', '+column)
                arrays[case['id']+'_'+column]=q[column]
            arrays[case['id']+'_t_s']=q['t_s']
        ax.set(xlabel=r'$t\ (\mathrm{ms})$',ylabel=ylabel);ax.legend(fontsize=7,frameon=False)
        save_figure(fig,figures/stem,arrays,data_output=data/(stem+'.npz'))
        _record(manifest,stem,'Unvalidated combined-arrangement prediction','Predicted quantities along integrated trajectories, evaluated with evolving in-plane gas-droplet slip. St_D_local uses tau_Stokes*norm(local mean gas velocity)/D0; the finite-Re coefficient timescale is distinct and exported in CSV. Stokes number alone does not determine transport or capture.')


def run(config_path=None):
    cfg_path=Path(config_path or PROJECT/'config/combined_flow_prediction.json');cfg=json.loads(cfg_path.read_text())
    if cfg.get('interaction_work_status')=='paused_for_plume_alignment_review':
        return dict(status='paused_for_plume_alignment_review',interaction_predictions_run=False)
    if cfg.get('co_directed',{}).get('enabled'):
        from .co_directed import run as run_upward
        return run_upward()
    out=(ROOT/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=True)
    if not out.is_relative_to(PROJECT.resolve()):raise ValueError('Output must stay in PZT analysis workspace')
    pzt,jet,xj,yj,D0,sources=sources_and_arrays(cfg);props=fluid_properties(cfg)
    authoritative=Path(cfg['pzt_statistics']);final=(ROOT/authoritative).parent.parent
    protected={str(f.resolve()):sha(f) for f in final.rglob('*') if f.is_file() and out not in f.resolve().parents and f.suffix.lower() in ('.npz','.csv','.png','.pdf','.svg')}
    for entry in sources['sources']:protected[entry['path']]=entry['sha256']
    manifest=[];peak=baseline(out,pzt,jet,xj,yj,ROOT/cfg['pzt_statistics'],manifest)
    # One configurable diameter population is shared by forces, Re and response/Stokes calculations.
    diameter=np.asarray(cfg['droplet_diameters_um'],float)
    if diameter.size:
        if not np.isfinite(diameter).all() or np.any(diameter<=0):raise ValueError('Positive finite configured diameters required')
        tau=relaxation_time(diameter*1e-6,props)
        write_csv(out/'data/droplet_model_parameters.csv',dict(configured_diameter_um=diameter,diameter_m=diameter*1e-6,
            tau_stokes_s=tau,tau_stokes_ms=tau*1e3))
    missing=missing_inputs(cfg);cases=[];case_records=[];convergence=[]
    if not missing:
        geometry=cfg['geometry'];position=np.asarray(geometry['jet_nozzle_position_in_pzt_frame_mm'],float)*1e-3
        angle=float(geometry['jet_rotation_degrees'])
        if abs(np.sin(np.deg2rad(angle)))>1e-10:raise ValueError('Proposed side-jet axis must be horizontal; document a different arrangement before using an inclined jet')
        if any(not np.isfinite(d) or d<=0 for d in cfg['droplet_diameters_um']):raise ValueError('Positive finite explicit diameters required')
        field=MeasuredGasField(xj,yj,jet['U'],jet['V'],position,angle)
        for r_index,release in enumerate(cfg['release_locations']):
            initial=release_state(release,pzt,ROOT/cfg['pzt_statistics'])
            for d_um in cfg['droplet_diameters_um']:
                ident=f'release_{r_index+1}_d{d_um:g}um';label=f'R{r_index+1}, d={d_um:g} um'
                try:field.sample(initial[:2])
                except UnsupportedPosition as error:
                    case_records.append(dict(id=ident,status='unsupported_release',reason=str(error),initial_state_SI=initial.tolist()));continue
                initial_re=float(reynolds(d_um*1e-6,initial[2:],field.sample(initial[:2]),props))
                if initial_re>=cfg['schiller_naumann_re_limit_exclusive']:
                    case_records.append(dict(id=ident,status='initial_Re_outside_drag_validity',initial_Re=initial_re,initial_state_SI=initial.tolist()));continue
                tr,status=integrate_droplet(field,initial,d_um*1e-6,props,cfg)
                fine,fine_status=integrate_droplet(field,initial,d_um*1e-6,props,cfg,step_scale=.5)
                # Compare at common supported time, independently of slightly different endpoint times.
                end=min(tr[-1,0],fine[-1,0]);common=np.linspace(0,end,201)
                errors=[float(np.max(np.abs(np.interp(common,tr[:,0],tr[:,k])-np.interp(common,fine[:,0],fine[:,k])))) for k in range(1,5)]
                if max(errors[:2])>cfg['solver']['position_convergence_tolerance_m'] or max(errors[2:])>cfg['solver']['velocity_convergence_tolerance_m_s']:
                    raise AssertionError('Prediction step convergence failed; refine configured time step')
                q=trajectory_quantities(field,tr,d_um*1e-6,props,D0,cfg)
                write_csv(out/'data'/(ident+'.csv'),q)
                cases.append(dict(id=ident,label=label,quantities=q))
                case_records.append(dict(id=ident,status=status,initial_state_SI=initial.tolist(),release=release,diameter_um=d_um,
                                         csv='data/'+ident+'.csv',end_time_s=float(tr[-1,0])))
                convergence.append(dict(id=ident,half_step_status=fine_status,errors_x_y_u_v=errors))
        if cases:prediction_figures(out,field,cases,manifest)
    status='awaiting_prediction_inputs' if missing else ('predictions_generated' if cases else 'no_supported_prediction_cases')
    receipt=dict(status=status,combined_experiment_conducted=False,missing_inputs=missing,
        configured_droplet_diameters_um=cfg['droplet_diameters_um'],diameter_source_or_assumption=cfg['diameter_source_or_assumption'],
        baseline_complete=True,prediction_case_count=len(cases),fixed_profile_maximum=peak,cases=case_records)
    write_json(out/'prediction_status.json',receipt);write_json(out/'configuration_used.json',cfg)
    write_json(out/'fluid_properties.json',props);write_json(out/'source_manifest.json',sources)
    write_json(out/'report_manifest.json',dict(figures=manifest,discovery='Explicit entries only; no recursive archive/old-theory discovery',
        text=['README.md','methods.tex','captions.tex'],archive_manifest=cfg['archive_manifest']))
    write_json(out/'validation/authoritative_file_hashes.json',protected)
    write_json(out/'validation/prediction_convergence.json',convergence)
    from .combined_flow_report import write_report
    write_report(out,cfg,sources,receipt,manifest)
    for path,digest in protected.items():assert sha(path)==digest,path
    from .prospective_studies import run as run_studies
    run_studies()  # Dimensionless calculations never depend on placement readiness.
    if cfg.get('design_study',{}).get('enabled'):
        receipt=json.loads((out/'prediction_status.json').read_text())
        write_report(out,cfg,sources,receipt,manifest)
    return receipt


def notebook_cells():
    import nbformat as nb
    active=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text())
    if active.get('interaction_work_status')=='paused_for_plume_alignment_review':
        from .plume_alignment import notebook_cells as alignment_cells
        return alignment_cells()
    if active.get('co_directed',{}).get('enabled'):
        from .co_directed import notebook_cells as upward_cells
        return upward_cells()
    text='''## Proposed side-jet combined-flow prediction

PZT-only measurements supply droplet velocities; the surrounding air was not
measured and the jet was absent. ZNMF-only measurements supply a separate gas
field. The combined experiment has not been conducted. Baselines are generated
without assuming alignment. The configured nominal diameter is a model estimate,
not a measurement: 5 micrometre aperture, assumed equal jet diameter, Rayleigh
factor 1.89 giving 9.45 micrometres, rounded to 10 micrometres for all calculations.
Diameter remains configurable. Predictions still require nozzle placement and
PZT release nodes in config/combined_flow_prediction.json.
The solver uses evolving vector slip in a prescribed measured mean gas field;
PZT velocities are used once at initialisation, never added to the gas field.
'''
    code='''from PZT_Atomiser_Analysis.analysis.combined_flow_prediction import run as run_prediction
prediction=run_prediction()
print(json.dumps(prediction,indent=2))
package=output/'Combined_Flow_Prediction'
manifest=json.loads((package/'report_manifest.json').read_text())
for figure in manifest['figures']:
    print(figure['kind'],figure['stem'])
    display(Image(filename=str(package/'figures'/(figure['stem']+'.png')),width=650))
'''
    cells=[nb.v4.new_markdown_cell(text),nb.v4.new_code_cell(code)]
    for c in cells:c.metadata['combined_flow_prediction']=True
    from .prospective_studies import notebook_cells as study_cells
    return cells+study_cells()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path)
    print(json.dumps(run(parser.parse_args().config),indent=2))
