"""Independent dimensionless audit from saved measurements and converged paths.

Run with --notebook to execute the audit and update only its notebook section.
No solver, calibration or placement changes are made by this module.
"""
import csv
import json
import os
import re
import shutil
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import numpy as np
from .combined_flow_prediction import ROOT, PROJECT, sha, write_json, write_csv
from .aligned_codirected import load_context, prepare
from .plotting import plt

RED='#d62728'; PURPLE='#7b2cbf'


def read(path):
    return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))


def rows_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    ctx=load_context();out=ctx['out'];dest=out/'Dimensionless_Numbers'
    legacy=['figures/'+stem+'.'+ext for stem in ('jet_local_Reynolds','jet_local_Strouhal','particle_Re','particle_Stokes') for ext in ('pdf','svg','png')]
    # Hash every existing scientific output except the exact generated files replaced.
    protected={str(p):sha(p) for p in (PROJECT/'Final_Results/Batch4_Cam1').rglob('*')
               if p.is_file() and 'Archive' not in p.parts and dest not in p.parents
               and p not in [out/x for x in legacy]+[out/n for n in ('dimensionless_manifest.json','report_manifest.json','README.md','captions.tex')]}
    if dest.exists() or any((out/x).exists() for x in legacy):
        archive=PROJECT/'Archive'/('Dimensionless_'+datetime.now(ZoneInfo('Australia/Sydney')).strftime('%Y%m%dT%H%M%S_%f_AEDT'))
        archive.mkdir(parents=True)
        records=[]
        for path in ([p for p in dest.rglob('*') if p.is_file()] if dest.exists() else [])+[out/x for x in legacy if (out/x).exists()]:
            rel=path.relative_to(out);target=archive/rel;target.parent.mkdir(parents=True,exist_ok=True)
            longpath=lambda p:('\\\\?\\'+str(p.resolve())) if os.name=='nt' else str(p)
            shutil.copy2(longpath(path),longpath(target));records.append(dict(original=str(rel),sha256=sha(path),operation='copy; original remains until regenerated'))
        write_json(archive/'archive_manifest.json',dict(scope='Only generated dimensionless outputs; measurements, source code and trajectories not archived',files=records))
    (dest/'figures').mkdir(parents=True,exist_ok=True);(dest/'data').mkdir(exist_ok=True)
    props=ctx['props'];rho=props['rho_air_kg_m3'];mu=props['mu_air_Pa_s'];rhol=props['rho_water_kg_m3'];nu=mu/rho
    D=float(ctx['gas_cfg']['diameter_mm'])*1e-3;f=float(ctx['gas_cfg']['forcing_Hz']);d=float(ctx['cfg']['droplet_diameters_um'][0])*1e-6
    assert D==.002 and f==120 and ctx['gas_cfg']['acquisition_Hz']==10
    tau=rhol*d*d/(18*mu);stkf=tau*f
    gas=ctx['gas'];X=gas['x_over_D0']*D*1e3;U=gas['Uc_m_per_s'].astype(float)
    valid=(gas['status']=='measured')&np.isfinite(X)&np.isfinite(U)&(U>0)
    Re=np.full(U.shape,np.nan);Sr=Re.copy();St=Re.copy()
    Re[valid]=U[valid]*D/nu;Sr[valid]=f*D/U[valid];St[valid]=tau*U[valid]/D
    np.testing.assert_allclose(X,ctx['sx'],atol=1e-12)
    np.testing.assert_allclose(St[valid]*Sr[valid],stkf,rtol=1e-14)
    # Independent exported-table check uses rho/mu rather than the nu formulation.
    previous=read(out/'data/jet_local_centreline_numbers.csv')
    np.testing.assert_allclose(previous['Re_D_c'],rho*D*np.where(valid,U,np.nan)/mu,equal_nan=True)
    np.testing.assert_allclose(previous['Sr_D_c'],Sr,equal_nan=True)
    np.testing.assert_allclose(previous['Stk_D_c'],St,equal_nan=True)
    write_csv(dest/'data/jet_local_centreline_numbers.csv',dict(native_index=np.arange(len(X)),x_over_D0=gas['x_over_D0'],distance_from_jet_nozzle_mm=X,Uc_m_s=U,source_status=gas['status'],defined_positive_measured=valid,Re_D_c=Re,Sr_D_c=Sr,Stk_D_c=St))
    scalar=dict(diameter_um=[d*1e6],tau_p_Stokes_s=[tau],tau_p_Stokes_ms=[tau*1e3],Stk_f=[stkf],D0_m=[D],forcing_Hz=[f],acquisition_Hz_NOT_USED=[ctx['gas_cfg']['acquisition_Hz']],rho_air_kg_m3=[rho],mu_air_Pa_s=[mu],nu_air_m2_s=[nu],rho_water_kg_m3=[rhol])
    write_csv(dest/'data/scalars.csv',scalar)
    field,states,initial,alignment,audit=prepare(ctx,write=False)
    summary=read(out/'data/trajectory_summary.csv');reference_number=int(summary['release_number'][summary['is_reference']][0]);paths=[];max_Re_error=0.;max_St_error=0.;gas_error=0.
    for row in summary:
        number=int(row['release_number']);q=read(out/'data'/f'release_{number:02d}.csv')
        gu=q['gas_u_p_m_s'];gv=q['gas_v_p_m_s'];su=gu-q['droplet_u_p_m_s'];sv=gv-q['droplet_v_p_m_s']
        slip=np.hypot(su,sv);gm=np.hypot(gu,gv);rp=rho*d*slip/mu;sp=tau*gm/D
        np.testing.assert_allclose(rp,q['Re_p'],rtol=1e-13,atol=1e-15);np.testing.assert_allclose(sp,q['Stk_D_local'],rtol=1e-13)
        np.testing.assert_allclose(q['ZNMF_downstream_mm'],q['h_mm']+alignment['entry']['ZNMF_downstream_mm'],atol=1e-12)
        # Resample the native two-component gas at every saved position, not Uc.
        sampled=np.array([field.sample(np.array([x,h])*1e-3) for x,h in zip(q['x_p_mm'],q['h_mm'])])
        np.testing.assert_allclose(sampled,np.c_[gu,gv],atol=2e-12)
        gas_error=max(gas_error,float(np.max(abs(sampled-np.c_[gu,gv]))))
        max_Re_error=max(max_Re_error,float(np.max(abs(rp-q['Re_p']))));max_St_error=max(max_St_error,float(np.max(abs(sp-q['Stk_D_local']))))
        F=np.where(rp<=ctx['cfg']['stokes_drag_re_threshold'],1,1+.15*rp**.687)
        columns=dict(t_s=q['t_s'],distance_from_plume_entry_h_mm=q['h_mm'],x_p_mm=q['x_p_mm'],distance_from_jet_nozzle_mm=q['ZNMF_downstream_mm'],gas_u_p_m_s=gu,gas_v_p_m_s=gv,droplet_u_p_m_s=q['droplet_u_p_m_s'],droplet_v_p_m_s=q['droplet_v_p_m_s'],slip_u_p_m_s=su,slip_v_p_m_s=sv,in_plane_slip_magnitude_m_s=slip,Re_p_in_plane=rp,Stk_D_path=sp,tau_p_Stokes_s=np.full(len(q),tau),tau_eff_coefficient_s=tau/F)
        write_csv(dest/'data'/f'path_{number:02d}_dimensionless.csv',columns)
        paths.append(dict(number=number,reference=number==reference_number,q=q,Re=rp,St=sp,gas_magnitude=gm))
    result=json.loads((out/'results_summary.json').read_text());reference=next(p for p in paths if p['reference']);q=reference['q'];minimum=result['reference_events']['first_velocity_minimum'];end=float(result['common_interval_h_mm'][1])
    events=[]
    for name,record in [('Entry',{k:float(q[k][0]) for k in ('t_s','h_mm','gas_u_p_m_s','gas_v_p_m_s','droplet_u_p_m_s','droplet_v_p_m_s')}),('Velocity minimum',minimum),('Comparison endpoint',{k:float(np.interp(end,q['h_mm'],q[k])) for k in ('t_s','h_mm','gas_u_p_m_s','gas_v_p_m_s','droplet_u_p_m_s','droplet_v_p_m_s')})]:
        # Endpoint components are interpolated on the saved path; event components
        # come from the supported Hermite root already exported by the solver.
        gu=record['gas_u_p_m_s'];gv=record['gas_v_p_m_s'];du=record['droplet_u_p_m_s'];dv=record['droplet_v_p_m_s'];rp=rho*d*np.hypot(gu-du,gv-dv)/mu
        events.append(dict(event=name,t_ms=record['t_s']*1e3,h_mm=record['h_mm'],ZNMF_distance_from_nozzle_mm=alignment['entry']['ZNMF_downstream_mm']+record['h_mm'],droplet_u_p_m_s=du,droplet_v_p_m_s=dv,gas_u_p_m_s=gu,gas_v_p_m_s=gv,slip_u_p_m_s=gu-du,slip_v_p_m_s=gv-dv,Re_p_in_plane=rp,Stk_D_path=tau*np.hypot(gu,gv)/D,Stk_f=stkf,tau_eff_coefficient_ms=tau*1e3/(1 if rp<=.01 else 1+.15*rp**.687)))
    rows_csv(dest/'data/reference_stations.csv',events)
    imax=int(np.nanargmax(np.where(valid,U,np.nan)));entryU=float(alignment['entry']['measured_Uc_m_s'])
    checks=[dict(station='Plume-entry centreline',Uc_m_s=entryU,distance_from_jet_nozzle_mm=alignment['entry']['ZNMF_downstream_mm'],Re_D_c=entryU*D/nu,Sr_D_c=f*D/entryU,Stk_D_c=tau*entryU/D),dict(station='Measured centreline maximum',Uc_m_s=float(U[imax]),distance_from_jet_nozzle_mm=float(X[imax]),Re_D_c=float(Re[imax]),Sr_D_c=float(Sr[imax]),Stk_D_c=float(St[imax]))]
    rows_csv(dest/'data/jet_reference_stations.csv',checks)
    assert abs(U[imax]-2.680478)<1e-5 and abs(gas['x_over_D0'][imax]-25.2096)<1e-6
    fig_stems=[];caps={}
    def save(fig,stem,caption):
        for ext in ('pdf','svg','png'):fig.savefig(dest/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig);fig_stems.append(stem);caps[stem]=caption
    for values,label,stem in [(Re,r'$Re_{D,c}$','jet_local_Reynolds'),(Sr,r'$Sr_{D,c}$','jet_local_Strouhal'),(St,r'$Stk_{D,c}$','jet_centreline_Stokes')]:
        fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained');ax.plot(X,values,color=RED,label='Measured centreline-based diagnostic')
        ax.set(xlabel=r'Distance from jet nozzle $X_j\ (\mathrm{mm})$',ylabel=label+' (–)');ax.legend(frameon=False,fontsize=8)
        save(fig,stem,'Local centreline-based quantity over the full valid measured ZNMF support; gaps and nonpositive Uc remain undefined. Distance is from the jet nozzle, not plume entry. D0=2 mm and forcing f=120 Hz; 10 Hz acquisition is excluded. Not a nozzle-exit quantity.')
    for key,label,stem in [('Re',r'$Re_p$','particle_Re'),('St',r'$Stk_{D,\mathrm{path}}$','particle_Stokes')]:
        fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained');first=True
        for p in paths:
            labelp='Reference droplet' if p['reference'] else ('Supported offset starts' if first else None)
            if not p['reference']:first=False
            ax.plot(p['q']['h_mm'],p[key],color=PURPLE if p['reference'] else '.65',lw=1.6 if p['reference'] else .8,label=labelp,zorder=4 if p['reference'] else 2)
        ax.axvline(end,color='.4',ls=':',lw=.8,label='End of PZT comparison support');ax.set(xlabel=r'Distance from plume entry $h\ (\mathrm{mm})$',ylabel=label+' (–)',xlim=(0,max(p['q']['h_mm'].max() for p in paths)));ax.legend(frameon=False,fontsize=8)
        save(fig,stem,'Every supported predicted path is shown to gas-support termination. Re_p uses both in-plane slip components at the droplet position; Stk_D,path uses the magnitude of local gas velocity and Stokes-limit tau_p. Off-axis gas is sampled from the full measured field, not substituted by Uc. Nominal 10 micrometre diameter is estimated, not measured. The vertical line marks the end of the shorter PZT-only comparison support.')
    fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained');ax.plot(q['h_mm'],tau*1e3/np.where(reference['Re']<=.01,1,1+.15*reference['Re']**.687),color=PURPLE,label=r'$\tau_{\mathrm{eff}}=\tau_p/F(Re_p)$');ax.axhline(tau*1e3,color='.4',ls='--',label=r'Stokes-limit $\tau_p$');ax.set(xlabel=r'Distance from plume entry $h\ (\mathrm{mm})$',ylabel=r'Response coefficient time $(\mathrm{ms})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'reference_response_times','tau_eff=tau_p/F is the instantaneous coefficient time in the existing finite-Re drag law, not a constant finite-Re relaxation eigen-time. F=1 at Re_p<=0.01, otherwise 1+0.15 Re_p^0.687. Stokes numbers use tau_p consistently.')
    # Retain established active report filenames while correcting their definitions.
    for rel in legacy:
        shutil.copy2(dest/rel,out/rel)
    def span(a):return f'{np.nanmin(a):.8g}–{np.nanmax(a):.8g}'
    allrp=np.concatenate([p['Re'] for p in paths]);allst=np.concatenate([p['St'] for p in paths])
    records=[
      dict(quantity='Re_D,c',formula='Uc D0 / nu_air',meaning='Local centreline jet inertia/viscosity',reference_scales='Uc(x), D0',value_or_range=span(Re),input_status='Measured sampled-mean Uc; documented D0; assumed air properties'),
      dict(quantity='Sr_D,c',formula='f D0 / Uc',meaning='Local centreline transit time / forcing period',reference_scales='D0/Uc(x), 1/f',value_or_range=span(Sr),input_status='Measured sampled-mean Uc; documented D0 and forcing f'),
      dict(quantity='tau_p (ms)',formula='rho_water d^2 / (18 mu_air)',meaning='Stokes-limit particle response time',reference_scales='Nominal spherical diameter d',value_or_range=f'{tau*1e3:.10g}',input_status='Model-estimated 10 um; recorded fluid-property assumptions'),
      dict(quantity='Stk_f',formula='tau_p f',meaning='Response time / forcing period',reference_scales='tau_p, 1/f',value_or_range=f'{stkf:.10g}',input_status='Estimated diameter; documented 120 Hz forcing'),
      dict(quantity='Stk_D,c',formula='tau_p Uc / D0',meaning='Response time / local centreline transit time',reference_scales='tau_p, D0/Uc(x)',value_or_range=span(St),input_status='Measured Uc; estimated diameter; recorded properties'),
      dict(quantity='Re_p',formula='rho_air d norm(u_g-v_d) / mu_air',meaning='In-plane slip-based particle Reynolds number',reference_scales='d, two-component vector slip',value_or_range=span(allrp),input_status='Measured prescribed gas; predicted droplet; estimated diameter'),
      dict(quantity='Stk_D,path',formula='tau_p norm(u_g(r_d)) / D0',meaning='Response time / local path gas transit time',reference_scales='tau_p, D0/norm(u_g)',value_or_range=span(allst),input_status='Measured field sampled along predicted paths; estimated diameter'),
      dict(quantity='Re0, St0',formula='U0 D0 / nu_air; f D0 / U0',meaning='Nozzle-reference quantities',reference_scales='Documented nozzle U0 required',value_or_range='Unresolved',input_status='No validated full-cycle nozzle-area mean square recovered'),
      dict(quantity='PZT-only measured Re_p',formula='Requires measured local gas minus droplet vector',meaning='Experimental slip Reynolds number',reference_scales='Surrounding gas velocity missing',value_or_range='Not identifiable',input_status='Droplet velocities alone do not measure gas slip')]
    rows_csv(dest/'data/summary.csv',records)
    evidence=PROJECT/'config/nozzle_reference_audit.json'
    nozzle=json.loads(evidence.read_text()) if evidence.exists() else dict(status='unresolved',missing='No documented phase-resolved nozzle-exit velocity over a forcing period or defensible expulsion-averaged reference speed in available experimental metadata.')
    write_json(dest/'nozzle_reference_audit.json',nozzle)
    table='| Quantity | Value/range |\n|---|---|\n'+'\n'.join(f"| {r['quantity']} | {r['value_or_range']} |" for r in records)
    text=f'''# Dimensionless-number audit

{table}

Ranges span the full valid measured jet centreline or all four complete supported paths, as applicable. No placement or trajectory was recalculated or changed. h=0 is the first supported measured plume-axis plane, placed at Xj={alignment['entry']['ZNMF_downstream_mm']:.9f} mm where interpolated Uc={entryU:.9f} m/s. The sampled centreline field is not rescaled to match that curve exactly.

At 20 C: rho_air={rho:.10g} kg/m3; mu_air={mu:.12g} Pa s; nu_air={nu:.12g} m2/s; rho_water={rhol:.12g} kg/m3. Atmospheric pressure 101325 Pa and dry-air properties are assumptions; humidity was not supplied. d={d*1e6:g} micrometres is a nominal model estimate: 1.89 times an assumed 5 micrometre emerging jet gives 9.45 micrometres, rounded to 10; neither jet/aperture equality nor measured droplet sizing is established.

The full valid centreline maximum is {U[imax]:.12g} m/s at x/D0={gas['x_over_D0'][imax]:.12g}, Xj={X[imax]:.12g} mm. Coordinates x/D0 are multiplied once by D0 in mm. Forcing is {f:g} Hz; acquisition 10 Hz is never used for Sr or Stk_f. Nonpositive velocities and source gaps remain undefined in local centreline diagnostics. No signed velocity is replaced by its magnitude. Norms occur only in vector-magnitude definitions of Re_p and Stk_D,path.

Stk_f uses forcing period 1/f; Stk_D,c uses centreline transit time D0/Uc; Stk_D,path uses local gas transit time D0/|u_g|. They answer different response questions. All use Stokes-limit tau_p. The solver retains finite-Re factor F; tau_eff=tau_p/F is a separately defined coefficient time, not an exact constant finite-Re relaxation time. Stk<1 does not prove entrainment, capture or successful delivery. The prescribed steady measured mean gas is one-way coupled; predicted paths are not measured combined-flow trajectories. PZT-only droplet velocities do not establish surrounding gas velocity or a measured slip Reynolds number.

Nozzle-reference finding: {nozzle['missing']} See nozzle_reference_audit.json for exact evidence and paper definition. Downstream maximum, entrance speed and signed full-cycle mean are not substituted for U0.

Reference-station values are in data/reference_stations.csv. Comparison endpoint uses component-wise interpolation on the saved path, followed by recalculation of norms; the refined minimum uses the existing supported event state. Path tables include every saved integration point through supported termination.

Reproduce without rerunning trajectories:
`.\\.venv\\Scripts\\python.exe -m PZT_Atomiser_Analysis.analysis.dimensionless_audit --notebook`
'''
    (dest/'README.md').write_text(text,encoding='utf-8')
    latex=r'''\section{Reynolds, Stokes and Strouhal numbers}
The local measured-centreline quantities are
\begin{equation}Re_{D,c}=\frac{U_c D_0}{\nu_g},\qquad Sr_{D,c}=\frac{fD_0}{U_c},\qquad U_c>0.\end{equation}
Source gaps are retained; these are not nozzle-exit quantities.
The nominal, model-estimated droplet diameter is $d=10\,\mu\mathrm m$;
the Rayleigh estimate assumes emerging jet diameter equals the $5\,\mu\mathrm m$ aperture and is not measured sizing.
\begin{equation}\tau_p=\frac{\rho_l d^2}{18\mu_g},\quad Stk_f=\tau_p f,\quad Stk_{D,c}=\frac{\tau_p U_c}{D_0},\quad Stk_{D,\mathrm{path}}=\frac{\tau_p|\mathbf u_g(\mathbf r_d)|}{D_0}.\end{equation}
The particle Reynolds estimate is in-plane and slip-based:
\begin{equation}Re_p=\frac{\rho_g d|\mathbf u_g(\mathbf r_d)-\mathbf v_d|}{\mu_g}.\end{equation}
Both gas components are interpolated within measured support at each predicted droplet position.
The measured steady mean gas is prescribed without two-way coupling. The existing finite-Re drag law is unchanged;
$\tau_{\mathrm{eff}}=\tau_p/F(Re_p)$ denotes only its local response-coefficient time.
The Stokes-limit response time defines all reported Stokes numbers, not a constant exact finite-Re relaxation.
The forcing period, centreline transit time and local path transit time are distinct reference scales;
$Stk<1$ alone proves neither complete entrainment nor successful delivery.
'''
    latex+=f'We use $D_0={D*1e3:g}$ mm and forcing $f={f:g}$ Hz, never the 10 Hz acquisition rate. At 20 degrees C, $\\rho_g={rho:.8g}$ kg/m$^3$, $\\mu_g={mu:.10g}$ Pa s, $\\nu_g={nu:.10g}$ m$^2$/s and $\\rho_l={rhol:.10g}$ kg/m$^3$. Pressure 101325 Pa and dry air are assumptions; humidity is unavailable. The resulting $\\tau_p={tau*1e3:.9f}$ ms and $Stk_f={stkf:.9f}$. '
    for r in checks:latex+=f"At $U_c={r['Uc_m_s']:.9f}$ m/s, $Re_{{D,c}}={r['Re_D_c']:.6f}$, $Sr_{{D,c}}={r['Sr_D_c']:.8f}$ and $Stk_{{D,c}}={r['Stk_D_c']:.8f}$. "
    latex+=f"The centreline maximum occurs at $X_j={X[imax]:.6f}$ mm ($x/D_0={gas['x_over_D0'][imax]:.6f}$). Across all supported predicted paths, $Re_p$ ranges from {allrp.min():.7g} to {allrp.max():.7g} and $Stk_{{D,\\mathrm{{path}}}}$ from {allst.min():.7g} to {allst.max():.7g}. "
    for r in events:latex+=f"{r['event']} at $h={r['h_mm']:.6f}$ mm has $Re_p={r['Re_p_in_plane']:.7f}$ and $Stk_{{D,\\mathrm{{path}}}}={r['Stk_D_path']:.7f}$. "
    latex+=r'''Cater and Soria (2002), equation (2.3), define the nozzle momentum-flow reference by
\begin{equation} U_0=\left[\frac{1}{A_0T_0}\int_0^{T_0}\int_{A_0}u_n^2\,\mathrm dA\,\mathrm dt\right]^{1/2},\qquad A_0=\pi D_0^2/4.\end{equation}
No representative full-cycle nozzle-area second moment was recovered, so $Re_0$ and $St_0$ are not assigned.
The unknown acquisition-phase mixture prevents assigning equal cycle weighting to the stored samples.
PZT-only droplet velocities do not establish ambient-gas slip.
'''
    latex+='\n\\begin{center}\\begin{tabular}{lrrr}\\hline\nReference station & $h$ (mm) & $Re_p$ & $Stk_{D,\\mathrm{path}}$ \\\\\n\\hline\n'
    for r in events:latex+=f"{r['event']} & {r['h_mm']:.6f} & {r['Re_p_in_plane']:.7f} & {r['Stk_D_path']:.7f} \\\\\n"
    latex+='\\hline\\end{tabular}\\end{center}\n'
    (dest/'methods_results.tex').write_text(latex,encoding='utf-8')
    (dest/'captions.tex').write_text('\n'.join('% '+k+'\n\\caption{'+v.replace('_',r'\_')+'}\n' for k,v in caps.items()),encoding='utf-8')
    for path,digest in protected.items():assert sha(path)==digest,path
    verification=dict(passed=True,protected_file_count=len(protected),saved_trajectories_and_Plume_Aligned_unchanged=True,calibration_not_reapplied=True,forcing_Hz=f,acquisition_excluded=True,coordinate_conversion='x_over_D0 * D0_m * 1000 exactly once',centreline_maximum=dict(Uc_m_s=float(U[imax]),x_over_D0=float(gas['x_over_D0'][imax]),Xj_mm=float(X[imax])),Stk_D_c_times_Sr_equals_Stk_f_max_error=float(np.max(abs(St[valid]*Sr[valid]-stkf))),maximum_path_Re_export_error=max_Re_error,maximum_path_Stk_export_error=max_St_error,maximum_resampled_gas_error_m_s=gas_error,diameter_scaling_verified=dict(Re_factor_if_doubled=2,tau_and_Stk_factor_if_doubled=4))
    np.testing.assert_allclose(rhol*(2*d)**2/(18*mu),4*tau)
    refslip=np.hypot(q['gas_u_p_m_s']-q['droplet_u_p_m_s'],q['gas_v_p_m_s']-q['droplet_v_p_m_s'])
    np.testing.assert_allclose(rho*(2*d)*refslip/mu,2*reference['Re'])
    exported=read(dest/'data/jet_local_centreline_numbers.csv')
    np.testing.assert_allclose(exported['Re_D_c'],Re,equal_nan=True);np.testing.assert_allclose(exported['Sr_D_c'],Sr,equal_nan=True)
    np.testing.assert_allclose(exported['Stk_D_c'],St,equal_nan=True)
    for p in paths:
        exported=read(dest/'data'/f'path_{p["number"]:02d}_dimensionless.csv')
        independent_Re=d*np.sqrt((exported['gas_u_p_m_s']-exported['droplet_u_p_m_s'])**2+(exported['gas_v_p_m_s']-exported['droplet_v_p_m_s'])**2)/nu
        independent_St=tau*np.sqrt(exported['gas_u_p_m_s']**2+exported['gas_v_p_m_s']**2)/D
        np.testing.assert_allclose(exported['Re_p_in_plane'],independent_Re,atol=1e-15)
        np.testing.assert_allclose(exported['Stk_D_path'],independent_St,atol=1e-15)
    verification['independent_exported_CSV_roundtrip_verified']=True
    write_json(dest/'verification.json',verification)
    write_json(dest/'source_manifest.json',dict(inputs=[dict(path=s['path'],sha256=s['sha256']) for s in ctx['sources']['sources']],trajectory_hashes={str(out/'data'/f'release_{int(r["release_number"]):02d}.csv'):sha(out/'data'/f'release_{int(r["release_number"]):02d}.csv') for r in summary},code=dict(path=str(Path(__file__)),sha256=sha(__file__)),nozzle_reference_evidence=dict(path=str(evidence),sha256=sha(evidence))))
    write_json(dest/'manifest.json',dict(figures=fig_stems,exports=['PDF','SVG','600 dpi PNG'],files=[str(p.relative_to(dest)) for p in dest.rglob('*') if p.is_file() and p.name!='manifest.json']))
    write_json(out/'dimensionless_manifest.json',dict(package='Dimensionless_Numbers',manifest='Dimensionless_Numbers/manifest.json',replaces_only=legacy))
    manifest=json.loads((out/'report_manifest.json').read_text());manifest['dimensionless_audit']=dict(package='Dimensionless_Numbers',manifest='Dimensionless_Numbers/manifest.json',LaTeX='Dimensionless_Numbers/methods_results.tex');write_json(out/'report_manifest.json',manifest)
    readme=out/'README.md';s=readme.read_text(encoding='utf-8').split('\n## Audited dimensionless results')[0]
    s+='\n## Audited dimensionless results\n\nSee [methods and summary](Dimensionless_Numbers/README.md), [complete definition table](Dimensionless_Numbers/data/summary.csv), [reference stations](Dimensionless_Numbers/data/reference_stations.csv), and [report LaTeX](Dimensionless_Numbers/methods_results.tex). Jet distances are from the nozzle; path distances are from plume entry. Full supported paths are shown, with the shorter PZT comparison endpoint marked. No experimental nozzle reference U0 was recovered.\n';readme.write_text(s,encoding='utf-8')
    cap=out/'captions.tex';s=cap.read_text(encoding='utf-8');new=(dest/'captions.tex').read_text(encoding='utf-8')
    for stem in ('jet_local_Reynolds','jet_local_Strouhal','particle_Re','particle_Stokes'):
        pattern=r'% '+re.escape(stem)+r'\n.*?(?=\n% |\Z)';replacement=re.search(pattern,new,re.S).group(0);s=re.sub(pattern,lambda m:replacement,s,flags=re.S)
    cap.write_text(s,encoding='utf-8')
    print(table);print(json.dumps(dict(jet_reference_stations=checks,reference_stations=events,verification=verification),indent=2))
    return dict(summary=records,jet_reference_stations=checks,reference_stations=events,figures=fig_stems)


def update_notebook():
    """Execute just the independent audit; preserve every pre-existing cell."""
    import asyncio, os, sys
    import nbformat as nb
    from nbclient import NotebookClient
    from jupyter_client import AsyncKernelManager
    from jupyter_client.kernelspec import KernelSpecManager
    target=PROJECT/'notebooks/PZT_final_analysis.ipynb';book=nb.read(target,as_version=4)
    cells=[nb.v4.new_markdown_cell('## Audited Reynolds, Stokes and Strouhal numbers\n\nLocal jet diagnostics use distance from the nozzle; path diagnostics use distance from plume entry. Nominal diameter is estimated. No nozzle-exit reference is invented.'),nb.v4.new_code_cell("from PZT_Atomiser_Analysis.analysis.dimensionless_audit import main as audit_dimensionless\ndimensionless_results=audit_dimensionless()\nfrom IPython.display import display,Image\nfrom pathlib import Path\nfor stem in dimensionless_results['figures']:\n    display(Image(filename=str(Path('PZT_Atomiser_Analysis/Final_Results/Batch4_Cam1/Combined_Flow_Prediction/Co_Directed/Dimensionless_Numbers/figures')/(stem+'.png')),width=700))")]
    for c in cells:c.metadata['dimensionless_audit']=True
    setup=nb.v4.new_code_cell(f'import sys\nsys.path.insert(0,{str(ROOT)!r})')
    mini=nb.v4.new_notebook(cells=[setup]+cells)
    os.environ['JUPYTER_RUNTIME_DIR']=str(PROJECT/'diagnostics/dimensionless_jupyter_runtime');Path(os.environ['JUPYTER_RUNTIME_DIR']).mkdir(exist_ok=True)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3');spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'];manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(mini,timeout=1800,resources={'metadata':{'path':str(ROOT)}});client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager);client.execute()
    original=[c for c in book.cells if not c.metadata.get('dimensionless_audit')]
    book.cells=original+mini.cells[1:];nb.write(book,target)
    for path in (PROJECT/'notebooks/Spatial_maps.ipynb',PROJECT/'Final_Results/Batch4_Cam1/methods/PZT_final_analysis.ipynb'):shutil.copy2(target,path)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--notebook',action='store_true');args=parser.parse_args()
    update_notebook() if args.notebook else main()
