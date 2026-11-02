"""Scientific outputs for the jet-under-plume entrance placement in (x_p,h)."""
import json
import numpy as np
from .plotting import plt
from .combined_flow_prediction import write_csv
from .aligned_codirected_report import clipped_curve
RED='#d62728';BLUE='#1f77b4';PURPLE='#7b2cbf'


def make_figures_and_report(ctx,field,histories,alignment,result):
    out=ctx['out'];lo,hi=result['common_interval_h_mm'];entry=alignment['entry'];p=ctx['profile'];sx=ctx['sx'];uc=ctx['uc'];good=ctx['valid'];s0=entry['ZNMF_downstream_mm'];mapped=sx-s0;reference=next(h for h in histories if h['reference']);q=reference['q'];n=reference['peak']+1;stems=[];captions={}
    def save(fig,stem,caption):
        for ext in ('pdf','svg','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig);stems.append(stem);captions[stem]=caption
    def axes(ax,label):ax.set(xlabel=r'$h\ (\mathrm{mm})$',ylabel=label,xlim=(lo,hi))
    first=True
    for h in histories:
        h['color']=PURPLE if h['reference'] else '.65';h['label']='Axis reference' if h['reference'] else (f'{len(histories)-1} supported offset starts' if first else None)
        if not h['reference']:first=False
    gx,gv=clipped_curve(mapped,uc,good,lo,hi);px,pv=clipped_curve(p['h_mm'],p['mean_v_p_m_s'],p['finite_pair'],lo,hi);rx,rv=clipped_curve(q['h_mm'][:n],q['droplet_v_p_m_s'][:n],np.ones(n,bool),lo,hi);_,rg=clipped_curve(q['h_mm'][:n],q['gas_v_p_m_s'][:n],np.ones(n,bool),lo,hi)
    fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained');ax.plot(gx,gv,color=RED,label='Measured ZNMF centreline; placed field');ax.plot(px,pv,color=BLUE,ls='--',label='Measured PZT-only fitted axis');ax.plot(0,entry['measured_Uc_m_s'],'ro',ms=3);axes(ax,r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'source_comparison_reproduced','Placement comparison: h=0 is the first valid measured PZT axis point, not a physical atomiser surface. The developing-branch jet station with interpolated centreline Uc=0.96 m/s is placed at h=0; Xj=x_j0+h. A crossing of these curves is not used to select the entrance.')
    for zoom,stem in [(False,'measured_and_predicted_axial_comparison'),(True,'reference_entry_response')]:
        fig,ax=plt.subplots(figsize=(7.4,4.8),layout='constrained');ax.plot(gx,gv,color=RED,label='Measured ZNMF centreline; rigid placement');ax.plot(px,pv,color=BLUE,ls='--',label='Measured PZT-only fitted axis');ax.plot(rx,rv,color=PURPLE,lw=1.7,label='Predicted reference droplet');ax.plot(rx,rg,color=RED,ls=':',label='Gas sampled on reference trajectory');axes(ax,r'Axial velocity $(\mathrm{m\,s^{-1}})$')
        if zoom:
            min_event=result['reference_events']['first_velocity_minimum'];end=min(hi,max(2.,2*min_event['h_mm'] if min_event else 2.));ax.set_xlim(0,end)
            visible=np.r_[gv[(gx>=0)&(gx<=end)],pv[(px>=0)&(px<=end)],rv[(rx>=0)&(rx<=end)],rg[(rx>=0)&(rx<=end)]]
            visible=visible[np.isfinite(visible)];span=max(float(np.ptp(visible)),.1);ax.set_ylim(float(visible.min())-.05*span,float(visible.max())+.08*span)
            for event,marker,label in [(result['reference_events']['first_axial_zero_slip'],'s','Zero axial slip'),(min_event,'o','Zero net acceleration / velocity minimum')]:
                if event:ax.plot(event['h_mm'],event['droplet_v_p_m_s'],marker,ms=4,mfc='white',mec='black',label=label)
        ax.legend(frameon=False,fontsize=7.5)
        save(fig,stem,'Reference droplet starts with measured near-base components, initially faster axially than the local gas. The later path is integrated from vector drag and effective gravity. The measured PZT-only ensemble profile is a comparison, never an imposed path. Zero axial slip is distinct from zero net acceleration and from zero vector slip. All main comparisons use the same supported h interval.')
    write_csv(out/'data/main_comparison_reference_prediction.csv',dict(h_mm=rx,predicted_v_p_m_s=rv,gas_sampled_v_p_m_s=rg));write_csv(out/'data/main_comparison_PZT_baseline.csv',dict(h_mm=px,mean_v_p_m_s=pv));write_csv(out/'data/main_comparison_gas_centreline.csv',dict(h_mm=gx,Uc_m_s=gv))
    xx,hh,ug,vg=field.common_mesh();limit=float(np.nanmax(abs(vg)))
    for detail,stem in [(False,'trajectories_upward_jet'),(True,'trajectories_entry_detail')]:
        fig,ax=plt.subplots(figsize=(6.4,6),layout='constrained');m=ax.contourf(xx*1e3,hh*1e3,np.ma.masked_invalid(vg),levels=np.linspace(-limit,limit,41),cmap='RdBu_r',corner_mask=False)
        for h in histories:
            v=h['q'];ax.plot(v['x_p_mm'],v['h_mm'],color=h['color'],label=h['label'],lw=1.7 if h['reference'] else .8,zorder=4 if h['reference'] else 2);ax.plot(v['x_p_mm'][0],0,'o',ms=3,color=h['color'])
        ax.axhline(0,color='.3',ls=':',lw=.8)
        if detail:
            top=min(hi,3.);xs=np.concatenate([h['q']['x_p_mm'][h['q']['h_mm']<=top] for h in histories]);ax.set(xlim=(xs.min()-.3,xs.max()+.3),ylim=(-.15,top+.15))
        ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$h\ (\mathrm{mm})$',aspect='equal');fig.colorbar(m,ax=ax,label=r'$\overline{v}_{g,p}\ (\mathrm{m\,s^{-1}})$',shrink=.8);fig.legend(*ax.get_legend_handles_labels(),loc='outside lower center',frameon=False,fontsize=8)
        save(fig,stem,'Jet-under-plume rigid placement. The entrance plane h=0 is indicated; the jet field extends below it. Only initially supported droplet starts are integrated. Genuine native gas masks and measured downstream centreline wandering are retained, with no row recentering, extrapolation or imposed outcome.')
    for key,label,stem in [('slip_v_p_m_s',r'$v_{g,p}-v_{d,p}\ (\mathrm{m\,s^{-1}})$','axial_slip'),('Re_p',r'$Re_p\ (-)$','particle_Re'),('Stk_D_local',r'$Stk_{D,\mathrm{local}}\ (-)$','particle_Stokes')]:
        fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
        for h in histories:
            v=h['q'];k=h['peak']+1;x,y=clipped_curve(v['h_mm'][:k],v[key][:k],np.ones(k,bool),lo,hi);ax.plot(x,y,color=h['color'],label=h['label'],lw=1.7 if h['reference'] else .8,zorder=4 if h['reference'] else 2)
        if key=='slip_v_p_m_s':ax.axhline(0,color='.5',lw=.6)
        axes(ax,label);ax.legend(frameon=False,fontsize=8)
        save(fig,stem,'Signed axial slip is gas minus droplet velocity. Re_p=rho_g*d*|u_g-v_d|/mu_g uses the full measured in-plane vector slip. Stk_D,local=tau_p*|u_g|/D0 uses the gas flow timescale D0/|u_g| and Stokes-limit tau_p=rho_l*d^2/(18*mu_g). Offset starts are unweighted sensitivity locations, not population samples.')
    common=np.genfromtxt(out/'data/reference_force_momentum_common_interval.csv',names=True,delimiter=',')
    for momentum,stem in [(False,'axial_acceleration'),(True,'axial_momentum_budget')]:
        fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
        fields=[('drag_impulse_y_p_N_s',RED,'-','Drag impulse'),('effective_gravity_impulse_y_p_N_s','.25','--','Effective-gravity impulse'),('momentum_change_y_p_kg_m_s',PURPLE,'-','Droplet momentum change')] if momentum else [('drag_a_y_p_m_s2',RED,'-','Aerodynamic drag'),('effective_gravity_a_y_p_m_s2','.25','--','Gravity minus buoyancy'),('net_a_y_p_m_s2',PURPLE,'-','Net droplet acceleration')]
        for key,color,style,label in fields:ax.plot(common['h_mm'],common[key]*(1e13 if momentum else 1),color=color,ls=style,label=label,lw=1.3)
        ax.axhline(0,color='.6',lw=.6);axes(ax,r'$\Delta p_h,J_h\ (10^{-13}\,\mathrm{kg\,m\,s^{-1}})$' if momentum else r'$a_h\ (\mathrm{m\,s^{-2}})$');ax.legend(frameon=False,fontsize=8)
        save(fig,stem,'Reference-droplet signed force/momentum balance from the measured base. Initially drag opposes upward motion because the droplet is faster than gas. Effective gravity is transformed into the plume frame. RK4-stage and independent saved-force trapezoid budgets are exported; no residual is set to zero. Momentum transfers through drag alone, not an additional jet force.')
    cut=np.genfromtxt(out/'data/release_coverage_profile.csv',names=True,delimiter=',');coverage=result['release_coverage'];fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained');ax.plot(cut['x_p_mm'],cut['mean_v_p_m_s'],color=BLUE,ls='--',label='Supported measured base-plane axial mean')
    for h in histories:ax.plot(h['q']['x_p_mm'][0],h['q']['droplet_v_p_m_s'][0],'o',color=h['color'],label=h['label'],ms=5)
    xmax=coverage['supported_native_plane_x_p_bounds_mm'][1];ax.axvspan(max(0,xmax),1.5,color='.93',hatch='//',label='No native PZT support on this plane');ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$',xlim=(-3,1.5));ax.legend(frameon=False,fontsize=7)
    save(fig,'initial_condition_coverage','The original seven offsets are audited on exactly h=0. Starts outside native PZT or gas support are excluded, not extrapolated. At this lowest plane the tilted native PZT window supports only one side of the fitted axis; no missing side, symmetry or complete plume width is invented. These are unweighted conditional starts.')
    props=ctx['props'];D0=ctx['D0'];tau=result['tau_stokes_s'];f=float(ctx['gas_cfg']['forcing_Hz']);nu=props['mu_air_Pa_s']/props['rho_air_kg_m3'];pos=good&(uc>0);v=np.where(pos,uc,np.nan)
    numbers=dict(ZNMF_downstream_mm=sx,h_mm=mapped,Uc_m_s=uc,status=ctx['gas']['status'],supported_positive_centreline=pos,Re_D_c=v*D0/nu,Sr_D_c=f*D0/v,Stk_D_c=tau*v/D0);write_csv(out/'data/jet_local_centreline_numbers.csv',numbers)
    write_csv(out/'data/scalar_quantities.csv',dict(nominal_estimated_diameter_um=[ctx['cfg']['droplet_diameters_um'][0]],tau_stokes_s=[tau],Stk_f=[tau*f],forcing_frequency_Hz=[f],acquisition_frequency_Hz_not_used=[ctx['gas_cfg']['acquisition_Hz']],D0_m=[D0],nu_air_m2_s=[nu],droplet_mass_kg=[props['rho_water_kg_m3']*np.pi*(ctx['cfg']['droplet_diameters_um'][0]*1e-6)**3/6]))
    for key,label,stem in [('Re_D_c',r'$Re_{D,c}\ (-)$','jet_local_Reynolds'),('Sr_D_c',r'$Sr_{D,c}\ (-)$','jet_local_Strouhal')]:
        fig,ax=plt.subplots(figsize=(6.5,4.4),layout='constrained');ax.plot(sx,numbers[key],color=RED);ax.set(xlabel=r'$X_j\ (\mathrm{mm})$',ylabel=label);save(fig,stem,'Measured local centreline definitions Re_D,c=Uc*D0/nu and Sr_D,c=f*D0/Uc retain missing gaps and exclude nonpositive Uc. These are not nozzle-exit quantities. Forcing frequency comes from source metadata, not the PIV acquisition rate.')
    result['dimensionless']=dict(tau_stokes_s=tau,Stk_f=tau*f,forcing_frequency_Hz=f,D0_m=D0,nu_air_m2_s=nu)
    text_report(ctx,result,captions)
    return stems


def text_report(ctx,result,captions):
    out=ctx['out'];e=result['entry'];r=result['reference_comparison'];events=result['reference_events'];minimum=events['minimum_over_common_interval'];zero=events['first_axial_zero_slip'];mom=result['momentum']['reference_common'];base=ctx['base'];candidate=json.loads((out/'interaction_plane_audit.json').read_text());rejected=candidate['rejected_release_count'];accepted=result['successfully_integrated_cases'];source=ctx['cfg']['diameter_source_or_assumption'];freq=ctx['gas_cfg']['forcing_Hz'];D0=ctx['D0']
    command=r'.\.venv\Scripts\python.exe -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook'
    text=f'''# Jet underneath the measured PZT plume: corrected base-entry prediction

Rebuild from the current PZT_droplets root:
```powershell
{command}
```
This regenerates the model and executed prediction notebook section, archives superseded generated Co_Directed outputs, and preserves source measurements, Plume_Aligned files and unrelated analyses. Current figures: report_manifest.json. All are title-free PDF, SVG and 600-dpi PNG; numerical traces and tables are in data/.

## Placement and input audit
This placement supersedes the mid-plume profile-crossing arrangement. No velocity-profile crossing selects the entrance. h=0 is the first valid saved fitted-axis PZT point: original laboratory (x,y)=({base['lab_x_mm']:.9f},{base['lab_y_mm']:.9f}) mm; saved plume y_p={base['y_p_mm']:.9f} mm. It is a measurement reference plane, not a documented physical atomiser surface. h=y_p-y_p_base. Saved velocities already use the 5.35 micrometres/pixel and 25 microsecond PZT calibration, which is not applied again and is not applied to ZNMF.

On the developing ZNMF branch, adjacent finite native samples bracket Uc=0.96 m/s. The interpolated station is x_j0={e['ZNMF_downstream_mm']:.9f} mm. Place its interpolated local centreline Yj={e['ZNMF_local_centreline_Y_mm']:.9f} mm on x_p=0 at h=0. One rigid placement gives x_p=Yc(x_j0)-Yj and h=Xj-x_j0, with velocity components u_g,p=-Vj and v_g,p=Uj. Thus Xj=x_j0+h everywhere; no coordinate stretch or downstream recentering is used. Original centreline motion and transverse structure remain. The jet-under-plume translation places its nozzle reference below h=0; no unmeasured nozzle or atomiser-surface travel is simulated. See entry_mapping.csv and transformation.json for both inverse maps.

The interpolated centreline curve is {e['measured_Uc_m_s']:.12f} m/s at h=0. The actual native bilinear field sample there is (u_g,p,v_g,p)=({e['sampled_gas_u_p_m_s']:.9f},{e['sampled_gas_v_p_m_s']:.9f}) m/s. Its axial difference from the curve is {e['sampled_minus_curve_v_m_s']:.9g} m/s because downstream interpolation and centreline extraction do not commute. This difference is explicitly retained; the field is neither rescaled nor made uniform. The 0.96 m/s condition applies to the entrance centreline only.

## Measured initial conditions and retained starts
The axis droplet starts with (u_d,p,v_d,p)=({e['PZT_axis_u_p_m_s']:.9f},{e['PZT_axis_v_p_m_s']:.9f}) m/s from the near-base aligned measurements. All starts use both local components from four finite native corners at the same h=0 plane. Of the original {len(candidate['release_candidates'])} offsets, {accepted} are supported and {rejected} are rejected. The tilted PZT support excludes the positive offsets at this lowest plane: their inverse-mapped laboratory heights are below the measured range. The plane is not moved upwards, the missing side is not mirrored or filled, and no extrapolation is used to retain seven starts. The retained starts are unweighted sensitivity locations, not a whole-plume distribution. Complete positions, source indices, interpolation weights and rejection reasons are exported in initial_conditions.csv, initial_interpolation_stencils.json and release_candidate_audit.csv. A full entry-plane FWHM cannot be established from the truncated support and is not reported as measured.

## Force law and scope
For constant spherical diameter, m=rho_l*pi*d^3/6, A=pi*d^2/4, w=u_g(r_d)-v_d and Re_p=rho_g*d*|w|/mu_g. Acceleration is F(Re_p)*w/tau_p+(1-rho_g/rho_l)*g_p, tau_p=rho_l*d^2/(18*mu_g). The existing documented hybrid is preserved: F=1 for Re_p<=0.01, otherwise F=1+0.15 Re_p^0.687, with Re_p<800 enforced at every RK stage. This stable response form equals the spherical C_D drag force without dividing by zero at zero slip. There is no additional jet force, velocity addition or imposed PZT-only deceleration.

Nominal diameter is {result['nominal_diameter_um']:g} micrometres, estimated rather than measured. Provenance: {source} Recorded fluid properties and sources are retained (20 C; assumed 101325 Pa; humidity unavailable). Constant diameter and no evaporation are assumptions. Gravity is {result['gravity_plume_m_s2']} m/s^2 in the inclined plume frame, not silently reset to a vertical remounting. The independent measured ZNMF mean is a steady, one-way-coupled prescribed gas: phase motion and gas modification by loading are unresolved. These are conditional calculated paths, not measured combined-flow or automatically ensemble-mean trajectories. The PZT-only ensemble profile is a comparison, not a later constraint or force, and does not establish stationary ambient air.

## Initial response, minimum and slip events
Initially v_d,p>v_g,p, so axial drag opposes upward droplet motion; effective gravity is also negative axially. The resolved reference minimum over the common measured interval is v_d,p={minimum['droplet_v_p_m_s']:.9f} m/s at h={minimum['h_mm']:.9f} mm, t={minimum['t_s']*1e3:.9f} ms. It is determined from net-acceleration sign changes with supported cubic-Hermite state interpolation and reevaluated gas/drag. The first negative-to-positive net-acceleration root is the onset of reacceleration; all roots and step-refinement differences are exported. An axial-slip root means v_g,p-v_d,p=0, not zero full vector slip. At a net-acceleration root, upward drag balances effective weight and axial slip is generally positive. reference_events.csv reports both slip components, vector magnitude and force terms at each event.
'''
    if zero:text+=f"First zero axial slip: h={zero['h_mm']:.9f} mm, t={zero['t_s']*1e3:.9f} ms, vector slip magnitude={zero['slip_magnitude_m_s']:.9g} m/s; net acceleration there={zero['net_a_y_p_m_s2']:.9f} m/s^2.\n"
    text+=f'''
## Outcomes over stated intervals
Over the common PZT/trajectory interval h={r['lower_h_mm']:.9f} to {r['upper_h_mm']:.9f} mm, the reference changes from {r['entry_droplet_v_p_m_s']:.9f} to {r['end_droplet_v_p_m_s']:.9f} m/s in {r['transit_time_ms']:.9f} ms, with lateral displacement {r['lateral_displacement_mm']:+.9f} mm. At the upper plane, sampled gas is {r['end_sampled_gas_v_p_m_s']:.9f} m/s and the measured PZT-only axis baseline is {r['PZT_only_axis_v_at_upper_m_s']:.9f} m/s. All supported trajectories also continue to their actual gas-support/time/validity termination, separately tabulated in trajectory_summary.csv and entry_exit_conditions.csv; these full-support endpoints are not confused with the end of the measured PZT comparison interval. Termination counts: {result['termination_reasons']}.

The reference common-interval axial momentum change is {mom['momentum_change_y_p_kg_m_s']:.12g} kg m/s, comprising drag impulse {mom['drag_impulse_y_p_N_s']:.12g} N s and effective-gravity impulse {mom['effective_gravity_impulse_y_p_N_s']:.12g} N s. RK4 residual {mom['momentum_residual_RK4_y_p_N_s']:.6g} N s; independent saved-force trapezoid residual {mom['momentum_residual_trapezoid_y_p_N_s']:.6g} N s. Stage closure is expected at roundoff and is not independent physical validation. No residual is forced to zero.

The full gas-supported reference endpoint is h={result['reference_full_support']['end_h_mm']:.9f} mm, v_d,p={result['reference_full_support']['end_v_p_m_s']:.9f} m/s, transit {result['reference_full_support']['end_time_s']*1e3:.9f} ms, and lateral displacement {result['reference_full_support']['delta_x_p_mm']:+.9f} mm. Its termination is {result['reference_full_support']['termination']}. This extends beyond the measured PZT-only comparison range.

Stk_D,local=tau_p*|u_g|/D0 uses flow time D0/|u_g|; Stk_f=tau_p*f. tau_p={result['tau_stokes_s']:.9g} s. D0={D0*1e3:g} mm and forcing f={freq:g} Hz come from ZNMF metadata, not its acquisition rate. Re_D,c=Uc*D0/nu and Sr_D,c=f*D0/Uc remain local centreline diagnostics, with missing/nonpositive Uc undefined. All gas, force and momentum arrays, response times, vector-slip Re and instantaneous balances are saved. Source hashes, transformations, support, analytical Stokes recovery, zero slip, drag direction, momentum closure and coarse/half-step convergence are checked in validation/.
'''
    (out/'README.md').write_text(text,encoding='utf-8')
    latex=r'''\section{Jet-under-plume entrance prediction}
The first valid measured fitted-axis point defines $h=0$, with $h=y_p-y_{p,0}$;
it is a measurement reference, not a documented atomiser surface.
A developing-branch station with measured centreline $U_c=0.96$ m/s is interpolated
only between adjacent finite samples. Its centreline is placed on the plume axis
at the entrance using one rigid map $h=X_j-x_{j0}$,
$x_p=Y_c(x_{j0})-Y_j$, and $(u_{g,p},v_{g,p})=(-V_j,U_j)$.
Both measured near-base droplet components initialise each supported release;
no profile crossing or mid-plume reset selects entry. Missing starts are rejected.
For $\mathbf w=\mathbf u_g(\mathbf r_d)-\mathbf v_d$,
\[
\dot{\mathbf r}_d=\mathbf v_d,\quad
\dot{\mathbf v}_d=F(Re_p)\mathbf w/\tau_p+(1-\rho_g/\rho_l)\mathbf g_p,
\quad Re_p=\rho_gd|\mathbf w|/\mu_g,\quad \tau_p=\rho_ld^2/(18\mu_g).
\]
This equals the spherical force with $m_d=\rho_l\pi d^3/6$ and $A_d=\pi d^2/4$.
The existing Stokes/Schiller--Naumann hybrid uses $F=1$ below $Re_p=0.01$ and
$F=1+0.15Re_p^{0.687}$ above, for $Re_p<800$; zero slip is handled safely.
The nominal $10\,\mu$m diameter is an estimate, not measured sizing.
Saved properties, individual experimental calibrations and rotated gravity are retained.
The independent measured gas mean is prescribed with steady one-way coupling;
phase motion, gas loading changes and evaporation are not resolved.
The PZT-only ensemble curve is a comparison, not a force or imposed trajectory.
'''
    latex+=f"The selected station is $x_{{j0}}={e['ZNMF_downstream_mm']:.6f}$ mm. Initial axial droplet and sampled-gas velocities are {e['PZT_axis_v_p_m_s']:.6f} and {e['sampled_gas_v_p_m_s']:.6f} m/s, so initial axial drag is negative. {accepted} of {len(candidate['release_candidates'])} candidate starts are supported. The reference minimum is {minimum['droplet_v_p_m_s']:.6f} m/s at $h={minimum['h_mm']:.6f}$ mm and $t={minimum['t_s']*1e3:.6f}$ ms. "
    if zero:latex+=f"Zero axial slip occurs at $h={zero['h_mm']:.6f}$ mm; it is distinct from zero net acceleration and does not imply zero vector slip. "
    latex+=f"At the common upper plane $h={r['upper_h_mm']:.6f}$ mm, reference velocity is {r['end_droplet_v_p_m_s']:.6f} m/s after {r['transit_time_ms']:.6f} ms, with transverse displacement {r['lateral_displacement_mm']:+.6f} mm. "
    full=result['reference_full_support']
    latex+=f"The full supported reference path ends at $h={full['end_h_mm']:.6f}$ mm with velocity {full['end_v_p_m_s']:.6f} m/s after {full['end_time_s']*1e3:.6f} ms; gas support, rather than the end of the PZT profile, terminates that path. "
    latex+=r'''Finite response produces slip; reacceleration begins when axial drag exceeds effective weight.
$Stk_{D,\rm local}=\tau_p|\mathbf u_g|/D_0$ uses time $D_0/|\mathbf u_g|$;
$Sr_{D,c}=fD_0/U_c$ is a distinct local jet Strouhal number.
Momentum satisfies $m_d\Delta\mathbf v_d=\int\mathbf F_Ddt+\int\mathbf F_{g,eff}dt+\boldsymbol\epsilon$,
checked by RK4-stage and independent saved-force quadrature and step refinement.
Calculated paths are not measured combined-flow results or automatically ensemble means;
no whole-plume distribution, capture efficiency or total plume momentum flux is inferred.
'''
    (out/'methods_results.tex').write_text(latex,encoding='utf-8')
    (out/'captions.tex').write_text('\n'.join('% '+stem+'\n\\caption{'+caption.replace('_',r'\_').replace('%',r'\%')+'}\n' for stem,caption in captions.items()),encoding='utf-8')


def render_saved():
    """Refresh final presentation from the saved converged base-entry trajectories."""
    from .aligned_codirected import load_context,prepare
    ctx=load_context();field,states,table,alignment,audit=prepare(ctx,write=False);out=ctx['out'];result=json.loads((out/'results_summary.json').read_text());histories=[]
    summary=np.genfromtxt(out/'data/trajectory_summary.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    for state in states:
        data=np.genfromtxt(out/'data'/f'release_{state["number"]:02d}.csv',names=True,delimiter=',');q={key:data[key] for key in data.dtype.names};summary_row=summary[summary['release_number']==state['number']][0]
        histories.append(dict(number=state['number'],reference=state['reference'],q=q,peak=int(np.argmax(q['h_mm'])),termination=str(summary_row['termination'])))
    if 'reference_full_support' not in result:
        row=summary[summary['is_reference']][0];result['reference_full_support']={key:(row[key].item() if hasattr(row[key],'item') else row[key]) for key in summary.dtype.names}
        (out/'results_summary.json').write_text(json.dumps(result,indent=2))
    make_figures_and_report(ctx,field,histories,alignment,result)
    return result
