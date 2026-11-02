"""Figures, source-coordinate audit and methods for aligned co-directed response."""
import json
import numpy as np
from .plotting import plt
from .combined_flow_prediction import write_csv,write_json
from .droplet_dynamics import relaxation_time
from .interaction_plane import supported_value

RED='#d62728';BLUE='#1f77b4';PURPLE='#7b2cbf'


def clipped_curve(x,y,valid,lower,upper):
    x=np.asarray(x);y=np.asarray(y);valid=np.asarray(valid)
    inside=(x>lower)&(x<upper);xx=np.r_[lower,x[inside],upper];yy=np.r_[np.nan,np.where(valid[inside],y[inside],np.nan),np.nan]
    for i,point in [(0,lower),(-1,upper)]:
        try:yy[i]=supported_value(x,y,valid,point)
        except ValueError:pass
    return xx,yy


def make_figures_and_report(ctx,field,histories,alignment,result):
    out=ctx['out'];p=ctx['profile'];cfg=ctx['cfg'];props=ctx['props'];D0=ctx['D0'];sx=ctx['sx'];uc=ctx['uc'];good=ctx['valid'];entry=alignment['entry'];cross=alignment['crossing'];lo,hi=result['common_interval_y_p_mm'];stems=[];captions={}
    def save(fig,stem,caption):
        for ext in ('pdf','svg','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig);stems.append(stem);captions[stem]=caption
    def axis(ax,ylabel):ax.set(xlabel=r'$y_p\ (\mathrm{mm})$',ylabel=ylabel,xlim=(lo,hi))
    colors=plt.cm.viridis(np.linspace(.05,.85,len(histories)))
    reference=next(h for h in histories if h['reference'])
    for h,color in zip(histories,colors):h['color']=PURPLE if h['reference'] else color;h['label']=rf'$x_{{p,0}}={h["q"]["x_p_mm"][0]:.3f}$ mm'
    paired_y=(sx+cfg['co_directed']['relative_axial_offset_mm']-ctx['origin'][1])/ctx['Q'][0,0]
    fig,ax=plt.subplots(figsize=(7.3,4.7),layout='constrained');domain=ctx['profile']['lab_y_mm'];comparison_domain=json.loads((out/'interaction_plane_audit.json').read_text())['original_comparison_domain_mm']
    gas_in=(sx>=comparison_domain[0])&(sx<=comparison_domain[1]);p_in=(p['lab_y_mm']>=comparison_domain[0])&(p['lab_y_mm']<=comparison_domain[1])
    ax.plot(paired_y[gas_in],np.where(good,uc,np.nan)[gas_in],color=RED,label='Measured ZNMF centreline; source pairing')
    ax.plot(p['y_p_mm'][p_in],np.where(p['finite_pair'],p['mean_v_p_m_s'],np.nan)[p_in],color=BLUE,ls='--',label='Measured PZT fitted axis')
    ax.plot(cross['y_p_mm'],cross['common_axial_velocity_m_s'],'ko',ms=3,label='Supported source-profile crossing');ax.axvline(entry['y_p_mm'],color='.35',ls=':',lw=.9,label='Selected native entry plane')
    ax.set(xlabel=r'$y_p\ (\mathrm{mm})$',ylabel=r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'source_comparison_reproduced','Corrected measured profiles retain the original pairing of ZNMF downstream distance with laboratory PZT height. Both source-pair coordinates are displayed in plume-axis units using the saved transform. The crossing is calculated only within adjacent finite support; the entry is a subsequent native gas station with actual sampled gas faster than the local PZT axis value. This source-pair mapping is not a spatial stretch of the gas used for integration.')
    q=reference['q'];n=reference['peak']+1;yy=q['y_p_mm'][:n];sample_x,sample_v=clipped_curve(yy,q['droplet_v_p_m_s'][:n],np.ones(n,bool),lo,hi);_,sample_g=clipped_curve(yy,q['gas_v_p_m_s'][:n],np.ones(n,bool),lo,hi)
    px,pv=clipped_curve(p['y_p_mm'],p['mean_v_p_m_s'],p['finite_pair'],lo,hi)
    placed_y=entry['y_p_mm']+sx-entry['ZNMF_downstream_mm'];gx,gv=clipped_curve(placed_y,uc,good,lo,hi)
    fig,ax=plt.subplots(figsize=(7.5,4.8),layout='constrained')
    ax.plot(gx,gv,color=RED,label='Measured ZNMF centreline; rigid placement');ax.plot(px,pv,color=BLUE,ls='--',label='Measured PZT-only fitted axis')
    ax.plot(sample_x,sample_v,color=PURPLE,lw=1.7,label='Predicted reference droplet');ax.plot(sample_x,sample_g,color=RED,ls=':',lw=1.5,label='Gas sampled along reference trajectory')
    axis(ax,r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'measured_and_predicted_axial_comparison','Reference prediction, measured PZT-only ensemble-axis baseline, measured gas centreline under one rigid placement, and actual gas sampled along the reference path. All are compared over the same supported ascending interval. The PZT ensemble profile is not an individual-droplet trajectory; the gas centreline is not necessarily the gas encountered by a laterally moving droplet.')
    write_csv(out/'data/main_comparison_reference_prediction.csv',dict(y_p_mm=sample_x,predicted_v_p_m_s=sample_v,gas_sampled_v_p_m_s=sample_g))
    write_csv(out/'data/main_comparison_PZT_baseline.csv',dict(y_p_mm=px,mean_v_p_m_s=pv));write_csv(out/'data/main_comparison_gas_centreline.csv',dict(y_p_mm=gx,Uc_m_s=gv))
    xx,yy,ug,vg=field.common_mesh();limit=float(np.nanmax(abs(vg)))
    for detail,stem in [(False,'trajectories_upward_jet'),(True,'trajectories_entry_detail')]:
        fig,ax=plt.subplots(figsize=(6.6,6),layout='constrained');m=ax.contourf(xx*1e3,yy*1e3,np.ma.masked_invalid(vg),levels=np.linspace(-limit,limit,41),cmap='RdBu_r',corner_mask=False)
        for h in histories:
            q=h['q'];ax.plot(q['x_p_mm'],q['y_p_mm'],color=h['color'],label=h['label'],lw=1.2);ax.plot(q['x_p_mm'][0],q['y_p_mm'][0],'o',ms=2.5,color=h['color'])
        ax.axhline(lo,color='.4',ls=':',lw=.7)
        if detail:
            allx=np.concatenate([h['q']['x_p_mm'][(h['q']['y_p_mm']>=lo)&(h['q']['y_p_mm']<=hi)] for h in histories]);ax.set(xlim=(allx.min()-.3,allx.max()+.3),ylim=(lo-.15,hi+.15))
        ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$y_p\ (\mathrm{mm})$',aspect='equal');fig.colorbar(m,ax=ax,label=r'$\overline{v}_{g,p}\ (\mathrm{m\,s^{-1}})$',shrink=.8);ax.legend(frameon=False,fontsize=7,loc='upper left')
        save(fig,stem,'Conditional trajectories begin with both locally measured PZT components on the supported entry plane. The measured two-component gas is rigidly placed once; measured jet wandering is retained. Contour masks are genuine measured gas support; integration stops at unsupported cells or domain boundaries. No atomiser-surface or near-nozzle release is modelled.')
    for key,label,stem,caption in [
        ('slip_v_p_m_s',r'$v_{g,p}-v_{d,p}\ (\mathrm{m\,s^{-1}})$','axial_slip','Signed gas-minus-droplet axial slip along the trajectories over their common supported ascending interval.'),
        ('Re_p',r'$Re_p\ (-)$','particle_Re','Particle Reynolds number uses the magnitude of the evolving two-component gas-minus-droplet slip, rho_air d |v_g-v_d| / mu_air. No out-of-plane motion is invented.'),
        ('Stk_D_local',r'$Stk_{D,\mathrm{local}}\ (-)$','particle_Stokes','Stk_D,local = tau_0 |v_g| / D0; tau_0 = rho_water d^2/(18 mu_air), and flow timescale D0/|v_g|. This is a Stokes-limit particle response diagnostic; finite-Re drag is used for the actual acceleration.'),
        ('net_a_y_p_m_s2',r'$a_{d,y_p}\ (\mathrm{m\,s^{-2}})$','axial_acceleration','Signed net axial acceleration includes finite-Re spherical drag and gravity-minus-buoyancy in the corrected inclined plume frame; no positive acceleration is imposed.')]:
        fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
        for h in histories:
            q=h['q'];n=h['peak']+1;xv,yv=clipped_curve(q['y_p_mm'][:n],q[key][:n],np.ones(n,bool),lo,hi);ax.plot(xv,yv,color=h['color'],label=h['label'],lw=1.2)
        if key in ('slip_v_p_m_s','net_a_y_p_m_s2'):ax.axhline(0,color='.5',lw=.6)
        axis(ax,label);ax.legend(frameon=False,fontsize=7,ncol=2);save(fig,stem,caption)
    tau=result['tau_stokes_s'];nu=props['mu_air_Pa_s']/props['rho_air_kg_m3'];f=float(ctx['gas_cfg']['forcing_Hz']);positive=good&(uc>0);uv=np.where(positive,uc,np.nan)
    numbers=dict(ZNMF_downstream_mm=sx,y_p_rigidly_placed_mm=placed_y,Uc_m_s=uc,status=ctx['gas']['status'],supported_positive_centreline=positive,Re_D_c=uv*D0/nu,Sr_D_c=f*D0/uv,Stk_D_c=tau*uv/D0)
    write_csv(out/'data/jet_local_centreline_numbers.csv',numbers)
    scalars=dict(nominal_estimated_diameter_um=[cfg['droplet_diameters_um'][0]],tau_stokes_s=[tau],Stk_f=[tau*f],forcing_frequency_Hz=[f],PIV_acquisition_frequency_Hz_not_used=[ctx['gas_cfg']['acquisition_Hz']],D0_m=[D0],nu_air_m2_s=[nu],rho_air_kg_m3=[props['rho_air_kg_m3']],rho_water_kg_m3=[props['rho_water_kg_m3']],mu_air_Pa_s=[props['mu_air_Pa_s']])
    write_csv(out/'data/scalar_quantities.csv',scalars)
    for key,label,stem in [('Re_D_c',r'$Re_{D,c}\ (-)$','jet_local_Reynolds'),('Sr_D_c',r'$Sr_{D,c}\ (-)$','jet_local_Strouhal')]:
        fig,ax=plt.subplots(figsize=(6.5,4.4),layout='constrained');ax.plot(sx,numbers[key],color=RED);ax.set(xlabel=r'$X_j\ (\mathrm{mm})$',ylabel=label);save(fig,stem,'Local measured-centreline quantity. Re_D,c = Uc D0/nu; Sr_D,c = f D0/Uc. D0 = 2 mm and measured jet actuation f = 120 Hz; PIV acquisition is 10 Hz and is not used as actuation. Gaps and nonpositive Uc remain undefined. These are not nozzle-exit Re0 or Sr0.')
    result['dimensionless']=dict(tau_stokes_s=tau,Stk_f=tau*f,D0_m=D0,forcing_frequency_Hz=f,acquisition_frequency_Hz_not_used=ctx['gas_cfg']['acquisition_Hz'],nu_air_m2_s=nu,definitions='Re_D,c=Uc D0/nu; Sr_D,c=f D0/Uc; Stk_D,c=tau0 Uc/D0; trajectory Stk_D,local=tau0 norm(gas)/D0; Stk_f=tau0 f')
    text_report(ctx,alignment,result,captions)
    return stems


def text_report(ctx,alignment,result,captions):
    out=ctx['out'];entry=result['entry'];cross=result['crossing'];r=result['reference_comparison'];gravity=result['gravity_plume_m_s2'];lo,hi=result['common_interval_y_p_mm'];numbers=result['dimensionless']
    command=r'.\.venv\Scripts\python.exe -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook'
    text=f'''# Co-directed response of the developed measured PZT plume

## Reproduce
From the local PZT_droplets workspace root:
```powershell
{command}
```
The command archives the old generated Co_Directed package, rebuilds the prediction, executes only its notebook cells, removes obsolete prediction displays and preserves unrelated notebook sections. Principal comparison: figures/measured_and_predicted_axial_comparison.pdf. Full figure list: report_manifest.json. Every figure is exported separately as PDF, SVG and 600-dpi PNG.

## Inputs and placement
Use the preserved full-support Plume_Aligned native grids and fitted-axis profile; report crops and the old fixed laboratory-column profile are not model inputs. Velocities are already SI, with 5.35 micrometres/pixel and 25 microseconds applied upstream; calibration is not applied again. The x_p/y_p frame follows the measured inclined plume, retaining the saved origin, orthonormal transform and inverse. Gravity in this frame is {gravity} m/s^2 (buoyancy reduces its magnitude by 1-rho_air/rho_water). Physical remounting is not asserted.

The source comparison originally paired ZNMF downstream Xj in mm with PZT laboratory height y_lab=Xj+offset, where offset={ctx['cfg']['co_directed']['relative_axial_offset_mm']} mm. That design pairing is retained to locate the interaction: the saved fitted-axis profile is indexed by its original lab_y_mm, and both comparison curves are displayed in y_p using the saved Q and reference. It is not evidence that experimental origins coincide. Changing the plume axial origin or report crop cannot relocate the paired laboratory height.

Supported gas-overtaking crossing: ZNMF Xj={cross['ZNMF_downstream_mm']:.9f} mm; PZT laboratory (x,y)=({cross['PZT_lab_x_mm']:.9f},{cross['PZT_lab_y_mm']:.9f}) mm; plume (x_p,y_p)=({cross['x_p_mm']:.9g},{cross['y_p_mm']:.9f}) mm; equal axial velocities {cross['common_axial_velocity_m_s']:.9f} m/s.

Selected entry is the first subsequent supported measured native gas station whose actual field sample exceeds the local axis droplet velocity: native column {entry['ZNMF_native_column']}, Xj={entry['ZNMF_downstream_mm']:.9f} mm; PZT laboratory ({entry['PZT_lab_x_mm']:.9f},{entry['PZT_lab_y_mm']:.9f}) mm; plume (0,{entry['y_p_mm']:.9f}) mm. Axis droplet v_p={entry['PZT_axis_v_p_m_s']:.9f} m/s; sampled gas v_p={entry['sampled_gas_v_p_m_s']:.9f} m/s. All candidates and valid brackets are in interaction_plane_audit.json.

At entry, align the measured jet centreline with x_p=0. Apply one rigid R90 placement of the jet local (Xj,Yj) coordinates and (Uj,Vj) components: u_g,p=-Vj, v_g,p=Uj. Downstream y_p=y_p,entry+Xj-Xj,entry and x_p=Yj,centreline(entry)-Yj. The gas centreline may wander thereafter; it is not recentered row-by-row. The source-pair comparison transformation has y_p=(y_lab-y_ref)/cos(theta) on the plume axis; this mapping is used to select entry only, not to stretch the gas field. The CSV ZNMF_profile_coordinate_mapping.csv explicitly lists both mappings. They agree at entry but are distinct downstream.

## Developed-droplet initial conditions
Representative transverse positions are the configured grid offsets times the saved rotated native horizontal-grid transverse spacing, all on the same fixed y_p entry plane. Both actual saved aligned mean components are interpolated using four finite original native PZT nodes; positions are inverse-mapped to the laboratory frame for interpolation. initial_conditions.csv records original coordinates, native indices, all four PZT weights, vector velocities, gas samples, slip and acceleration; initial_interpolation_stencils.json includes PZT and gas source corner values and weights. These conditional predictions start with already-developed droplets, not at the atomiser or nozzle. Mean measured components initialise deterministic reference trajectories; the ensemble baseline is not a measured individual trajectory.

## Dynamics and assumptions
Nominal constant diameter is 10 micrometres, a model estimate motivated by a 5-micrometre aperture and classical ratio 1.89 (9.45 micrometres rounded), assuming the emerging jet diameter equals the aperture. Neither diameter nor size distribution was measured. Recorded water/air properties at 20 C and assumed 101325 Pa are retained; humidity is unavailable and evaporation is neglected. Source references and values remain in source_manifest.json and the preserved fluid_properties.json.

The independently measured two-component ZNMF mean is prescribed as the gas. One-way droplet momentum coupling, a two-dimensional deterministic mean-field trajectory and no turbulence/evaporation/size distribution are model assumptions. Velocities are never added, and downstream PZT-only velocities are not forces or imposed trajectory values. Acceleration is (F(Re_p)/tau0)(v_g-v_d)+(1-rho_air/rho_water)g_p; tau0=rho_water*d^2/(18*mu_air). Existing spherical drag uses exact Stokes F=1 for Re_p<=0.01, otherwise Schiller-Naumann F=1+0.15 Re_p^0.687 for Re_p<800. Re_p uses the two-component slip magnitude. Integration stops at unsupported four-corner gas cells or measured-domain boundaries, with no extrapolation. No positive velocity, acceleration or transport outcome is enforced.

## Common-interval result
Successfully integrated cases: {result['successfully_integrated_cases']}. All releases and baselines are compared over the same supported ascending plume interval [{lo:.9f},{hi:.9f}] mm. Reference droplet: entry v_p={r['entry_droplet_v_p_m_s']:.9f} m/s; final v_p={r['end_droplet_v_p_m_s']:.9f} m/s; change {r['delta_v_p_m_s']:+.9f} m/s; transit {r['transit_time_ms']:.6f} ms; transverse displacement {r['lateral_displacement_mm']:+.6f} mm. At the same upper plane, sampled gas v_p={r['end_sampled_gas_v_p_m_s']:.9f} m/s and the measured PZT-only ensemble axis baseline is {r['PZT_only_axis_v_at_upper_m_s']:.9f} m/s. trajectory_summary.csv records full-path termination and common_height_comparison.csv records comparable interval results. Termination counts: {result['termination_reasons']}.

## Dimensionless response
Stokes-limit response tau0={numbers['tau_stokes_s']:.9g} s. Trajectory Stk_D,local=tau0*|v_g|/D0 uses flow time D0/|v_g|, not slip or forcing period; finite-Re coefficient response tau0/F(Re_p) is also exported. Stk_f=tau0*f={numbers['Stk_f']:.9g}. Local centreline diagnostics Re_D,c=Uc*D0/nu, Sr_D,c=f*D0/Uc and Stk_D,c=tau0*Uc/D0 use measured positive Uc only; missing stations remain gaps. D0=2 mm, actuation f=120 Hz from ZNMF configuration, not the 10-Hz PIV acquisition rate. These are local centreline-based quantities, not nozzle-exit Re0/Sr0. Stk is Stokes number and Sr is Strouhal number.

## Verification and limitations
Source and all Plume_Aligned files are hash-protected. Position/vector inverse transforms, unchanged sample counts and masks, gravity direction, finite interpolation stencils, actual gas support, slip-direction drag, and coarse/half-step convergence are checked. validation/convergence.json records errors and termination. Discrete starting locations and deterministic mean velocities do not imply measured particle histories, collection efficiency or an optimal arrangement. Out-of-plane motion, changing diameters, unresolved gas fluctuations and two-way coupling are not measured/modelled. Gas sampled along a drifting trajectory differs from the measured centreline; the two are separately plotted and exported.
'''
    (out/'README.md').write_text(text,encoding='utf-8')
    # Standalone concise LaTeX section; all detailed audit tables remain linked in README.
    latex=r'''\section{Co-directed developed-droplet response}
The preserved plume-aligned PZT means initialise both droplet components at a supported entry plane.
The original design pairing $X_j=y_{\rm lab}-\delta$ is retained when finding the crossing;
$y_p=(y_{\rm lab}-y_{\rm ref})/\cos\theta$ on the fitted axis. It does not imply coincident experimental origins.
One rigid placement then aligns the measured jet centreline with the fitted plume axis at entry;
there is no downstream recentering or rescaling. The independent measured two-component gas mean
is prescribed, without velocity addition, two-way coupling or imposed downstream PZT velocities.
Nominal $d=10\,\mu$m is an aperture/Rayleigh-based estimate, not measured sizing.
The force law is
\[\dot{\mathbf v}_d=F(Re_p)(\mathbf v_g-\mathbf v_d)/\tau_0+(1-\rho_a/\rho_w)\mathbf g_p,\qquad
\tau_0=\rho_w d^2/(18\mu_a).\]
Existing Stokes/Schiller--Naumann drag is used with the saved fluid properties, constant diameter
and no evaporation. The saved gravity vector is transformed into the measured inclined plume frame.
Four-finite-corner interpolation never extrapolates; trajectories stop at missing gas support.
The measured PZT-only ensemble axis curve is a baseline, not an individual-droplet trajectory.
'''
    latex+=f'Crossing: $X_j={cross["ZNMF_downstream_mm"]:.6f}$ mm; entry: $X_j={entry["ZNMF_downstream_mm"]:.6f}$ mm, $y_p={entry["y_p_mm"]:.6f}$ mm. '
    latex+=f'Over $y_p={lo:.6f}$--${hi:.6f}$ mm, the reference droplet changes from {r["entry_droplet_v_p_m_s"]:.6f} to {r["end_droplet_v_p_m_s"]:.6f} m/s in {r["transit_time_ms"]:.6f} ms, with transverse displacement {r["lateral_displacement_mm"]:+.6f} mm. '
    latex+=r'''Particle $Re_p=\rho_a d|\mathbf v_g-\mathbf v_d|/\mu_a$ uses vector slip;
$Stk_{D,\rm local}=\tau_0|\mathbf v_g|/D_0$ uses flow time $D_0/|\mathbf v_g|$.
Local centreline $Re_{D,c}=U_cD_0/\nu_a$ and $Sr_{D,c}=fD_0/U_c$ retain gaps and exclude nonpositive $U_c$;
$D_0=2$ mm and actuation $f=120$ Hz, not the 10-Hz acquisition rate. These are not nozzle-exit quantities.
No collection efficiency, experimentally validated trajectory or optimal arrangement is claimed.
'''
    (out/'methods_results.tex').write_text(latex,encoding='utf-8')
    def escape(text):return text.replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
    (out/'captions.tex').write_text('\n'.join('% '+stem+'\n\\caption{'+escape(caption)+'}\n' for stem,caption in captions.items()),encoding='utf-8')
