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
    first_offset=True
    for h in histories:
        h['color']=PURPLE if h['reference'] else '0.65'
        h['label']='Axis reference' if h['reference'] else ('Six unweighted sensitivity starts' if first_offset else None)
        if not h['reference']:first_offset=False
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
            q=h['q'];ax.plot(q['x_p_mm'],q['y_p_mm'],color=h['color'],label=h['label'],lw=1.8 if h['reference'] else .8,zorder=4 if h['reference'] else 2);ax.plot(q['x_p_mm'][0],q['y_p_mm'][0],'o',ms=2.5,color=h['color'])
        ax.axhline(lo,color='.4',ls=':',lw=.7)
        if detail:
            detail_upper=min(lo+3.,hi)
            allx=np.concatenate([h['q']['x_p_mm'][(h['q']['y_p_mm']>=lo)&(h['q']['y_p_mm']<=detail_upper)] for h in histories]);ax.set(xlim=(allx.min()-.3,allx.max()+.3),ylim=(lo-.15,detail_upper+.15))
        ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$y_p\ (\mathrm{mm})$',aspect='equal')
        from matplotlib.ticker import MaxNLocator
        cb=fig.colorbar(m,ax=ax,label=r'$\overline{v}_{g,p}\ (\mathrm{m\,s^{-1}})$',shrink=.8);cb.locator=MaxNLocator(nbins=6);cb.update_ticks()
        fig.legend(*ax.get_legend_handles_labels(),frameon=False,fontsize=8,ncol=1,loc='outside lower center')
        save(fig,stem,('Entry detail shows the first 3 mm of axial development. ' if detail else '')+'Conditional trajectories begin with both locally measured PZT components on the supported entry plane. The measured two-component gas is rigidly placed once; measured jet wandering is retained. Contour masks are genuine measured gas support; integration stops at unsupported cells or domain boundaries. No atomiser-surface or near-nozzle release is modelled.')
    for key,label,stem,caption in [
        ('slip_v_p_m_s',r'$v_{g,p}-v_{d,p}\ (\mathrm{m\,s^{-1}})$','axial_slip','Signed gas-minus-droplet axial slip along the trajectories over their common supported ascending interval.'),
        ('Re_p',r'$Re_p\ (-)$','particle_Re','Particle Reynolds number uses the magnitude of the evolving two-component gas-minus-droplet slip, rho_air d |v_g-v_d| / mu_air. No out-of-plane motion is invented.'),
        ('Stk_D_local',r'$Stk_{D,\mathrm{local}}\ (-)$','particle_Stokes','Stk_D,local = tau_0 |v_g| / D0; tau_0 = rho_water d^2/(18 mu_air), and flow timescale D0/|v_g|. This is a Stokes-limit particle response diagnostic; finite-Re drag is used for the actual acceleration.'),
        ]:
        fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
        for h in histories:
            q=h['q'];n=h['peak']+1;xv,yv=clipped_curve(q['y_p_mm'][:n],q[key][:n],np.ones(n,bool),lo,hi);ax.plot(xv,yv,color=h['color'],label=h['label'],lw=1.8 if h['reference'] else .8,zorder=4 if h['reference'] else 2)
        if key in ('slip_v_p_m_s','net_a_y_p_m_s2'):ax.axhline(0,color='.5',lw=.6)
        axis(ax,label);ax.legend(frameon=False,fontsize=7,ncol=2);save(fig,stem,caption)
    # Reference force balance, with signed gravity, not an extra jet force.
    common=np.genfromtxt(out/'data/reference_force_momentum_common_interval.csv',names=True,delimiter=',')
    fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
    for key,color,style,label in [('drag_a_y_p_m_s2',RED,'-','Aerodynamic drag'),('effective_gravity_a_y_p_m_s2','0.25','--','Gravity minus buoyancy'),('net_a_y_p_m_s2',PURPLE,'-','Net droplet acceleration')]:
        ax.plot(common['y_p_mm'],common[key],color=color,ls=style,label=label,lw=1.4)
    ax.axhline(0,color='.6',lw=.6);axis(ax,r'$a_{y_p}\ (\mathrm{m\,s^{-2}})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'axial_acceleration','Reference-droplet force balance: signed aerodynamic drag acceleration, effective gravity and their sum. Upward drag need not imply net upward acceleration: it must exceed effective weight. Gas momentum transfer is already the aerodynamic drag term; there is no additional jet force.')
    fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
    for key,color,style,label in [('drag_impulse_y_p_N_s',RED,'-','Drag impulse'),('effective_gravity_impulse_y_p_N_s','0.25','--','Effective-gravity impulse'),('momentum_change_y_p_kg_m_s',PURPLE,'-','Droplet momentum change')]:
        ax.plot(common['y_p_mm'],common[key]*1e13,color=color,ls=style,label=label,lw=1.4)
    ax.axhline(0,color='.6',lw=.6);axis(ax,r'$\Delta p_{y_p},J_{y_p}\ (10^{-13}\,\mathrm{kg\,m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
    save(fig,'axial_momentum_budget','Single model-droplet axial momentum change equals integrated drag plus integrated effective gravity. RK4-stage impulses and residuals are exported, together with an independent trapezoidal check from saved forces; neither residual is forced to zero. These are not whole-plume momentum fluxes.')
    coverage=result['release_coverage'];cut=np.genfromtxt(out/'data/release_coverage_profile.csv',names=True,delimiter=',')
    fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained');ax.plot(cut['x_p_mm'],cut['mean_v_p_m_s'],color=BLUE,ls='--',label='Measured entry-plane axial mean')
    for h in histories:
        q=h['q'];ax.plot(q['x_p_mm'][0],q['droplet_v_p_m_s'][0],'o',color=h['color'],ms=5,label=h['label'])
    ax.plot([coverage['left_halfmax_crossing_x_p_mm'],coverage['right_halfmax_crossing_x_p_mm']],[coverage['halfmax_threshold_v_p_m_s']]*2,'k:',lw=1,label='Measured half-maximum width')
    extent=np.ceil(max(abs(coverage['left_halfmax_crossing_x_p_mm']),abs(coverage['right_halfmax_crossing_x_p_mm']))*1.6*2)/2
    ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$',xlim=(-extent,extent));ax.legend(frameon=False,fontsize=8)
    save(fig,'initial_condition_coverage','Audit of the retained seven initial conditions: one axis reference and six unweighted central-region sensitivity starts. Their regular spacing was a numerical convenience, not a measured particle distribution. The measured entry-plane profile and its asymmetric half-maximum crossings show their limited coverage. No population weights, whole-plume representativeness or collection efficiency are inferred.')
    tau=result['tau_stokes_s'];nu=props['mu_air_Pa_s']/props['rho_air_kg_m3'];f=float(ctx['gas_cfg']['forcing_Hz']);positive=good&(uc>0);uv=np.where(positive,uc,np.nan)
    numbers=dict(ZNMF_downstream_mm=sx,y_p_rigidly_placed_mm=placed_y,Uc_m_s=uc,status=ctx['gas']['status'],supported_positive_centreline=positive,Re_D_c=uv*D0/nu,Sr_D_c=f*D0/uv,Stk_D_c=tau*uv/D0)
    write_csv(out/'data/jet_local_centreline_numbers.csv',numbers)
    scalars=dict(nominal_estimated_diameter_um=[cfg['droplet_diameters_um'][0]],tau_stokes_s=[tau],Stk_f=[tau*f],forcing_frequency_Hz=[f],PIV_acquisition_frequency_Hz_not_used=[ctx['gas_cfg']['acquisition_Hz']],D0_m=[D0],nu_air_m2_s=[nu],rho_air_kg_m3=[props['rho_air_kg_m3']],rho_water_kg_m3=[props['rho_water_kg_m3']],mu_air_Pa_s=[props['mu_air_Pa_s']])
    scalars.update(droplet_mass_kg=[props['rho_water_kg_m3']*np.pi*(cfg['droplet_diameters_um'][0]*1e-6)**3/6],projected_area_m2=[np.pi*(cfg['droplet_diameters_um'][0]*1e-6)**2/4])
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
Existing unweighted numerical sensitivity positions are the configured grid offsets times the saved rotated native horizontal-grid transverse spacing, all on the same fixed y_p entry plane. Both actual saved aligned mean components are interpolated using four finite original native PZT nodes; positions are inverse-mapped to the laboratory frame for interpolation. initial_conditions.csv records original coordinates, native indices, all four PZT weights, vector velocities, gas samples, slip and acceleration; initial_interpolation_stencils.json includes PZT and gas source corner values and weights. At the axis entry, interpolation of the saved one-dimensional axis profile gives {entry['PZT_saved_axis_profile_interpolated_v_p_m_s']:.9f} m/s, versus {entry['PZT_axis_v_p_m_s']:.9f} m/s from the actual native-cell initial-condition interpolation (difference {entry['PZT_axis_native_field_minus_profile_m_s']:.9g} m/s). Sampling a bilinear field along an inclined line and then interpolating the resulting one-dimensional samples is not identical to evaluating the native bilinear cell at the final point; no velocity matching is forced. These conditional predictions start with already-developed droplets, not at the atomiser or nozzle. Mean measured components initialise deterministic reference trajectories; the ensemble baseline is not a measured individual trajectory.

## Dynamics and assumptions
Nominal constant diameter is 10 micrometres, a model estimate motivated by a 5-micrometre aperture and classical ratio 1.89 (9.45 micrometres rounded), assuming the emerging jet diameter equals the aperture. Neither diameter nor size distribution was measured. Recorded water/air properties at 20 C and assumed 101325 Pa are retained; humidity is unavailable and evaporation is neglected. Source references and values remain in source_manifest.json and the preserved fluid_properties.json.

The independently measured two-component ZNMF mean is prescribed as a steady gas field. This is a steady, one-way-coupled mean-field approximation: droplets respond to the gas, gas modification by droplet loading is not calculated, and phase-dependent ZNMF motion is unresolved. A two-dimensional deterministic mean-field trajectory and no turbulence/evaporation/size distribution are model assumptions. Velocities are never added, and downstream PZT-only velocities are not forces or imposed trajectory values. Acceleration is (F(Re_p)/tau0)(v_g-v_d)+(1-rho_air/rho_water)g_p; tau0=rho_water*d^2/(18*mu_air). Existing spherical drag uses exact Stokes F=1 for Re_p<=0.01, otherwise Schiller-Naumann F=1+0.15 Re_p^0.687 for Re_p<800. Re_p uses the two-component slip magnitude. Integration stops at unsupported four-corner gas cells or measured-domain boundaries, with no extrapolation. No positive velocity, acceleration or transport outcome is enforced.

## Common-interval result
Successfully integrated cases: {result['successfully_integrated_cases']}. All releases and baselines are compared over the same supported ascending plume interval [{lo:.9f},{hi:.9f}] mm. Reference droplet: entry v_p={r['entry_droplet_v_p_m_s']:.9f} m/s; final v_p={r['end_droplet_v_p_m_s']:.9f} m/s; change {r['delta_v_p_m_s']:+.9f} m/s; transit {r['transit_time_ms']:.6f} ms; transverse displacement {r['lateral_displacement_mm']:+.6f} mm. At the same upper plane, sampled gas v_p={r['end_sampled_gas_v_p_m_s']:.9f} m/s and the measured PZT-only ensemble axis baseline is {r['PZT_only_axis_v_at_upper_m_s']:.9f} m/s. trajectory_summary.csv records full-path termination and common_height_comparison.csv records comparable interval results. Termination counts: {result['termination_reasons']}.

## Dimensionless response
Stokes-limit response tau0={numbers['tau_stokes_s']:.9g} s. Trajectory Stk_D,local=tau0*|v_g|/D0 uses flow time D0/|v_g|, not slip or forcing period; finite-Re coefficient response tau0/F(Re_p) is also exported. Stk_f=tau0*f={numbers['Stk_f']:.9g}. Local centreline diagnostics Re_D,c=Uc*D0/nu, Sr_D,c=f*D0/Uc and Stk_D,c=tau0*Uc/D0 use measured positive Uc only; missing stations remain gaps. D0=2 mm, actuation f=120 Hz from ZNMF configuration, not the 10-Hz PIV acquisition rate. These are local centreline-based quantities, not nozzle-exit Re0/Sr0. Stk is Stokes number and Sr is Strouhal number.

## Verification and limitations
Source and all Plume_Aligned files are hash-protected. Position/vector inverse transforms, unchanged sample counts and masks, gravity direction, finite interpolation stencils, actual gas support, slip-direction drag, and coarse/half-step convergence are checked. validation/convergence.json records errors and termination. Discrete starting locations and deterministic mean velocities do not imply measured particle histories, collection efficiency or an optimal arrangement. Out-of-plane motion, changing diameters, unresolved gas fluctuations and two-way coupling are not measured/modelled. Gas sampled along a drifting trajectory differs from the measured centreline; the two are separately plotted and exported.
'''
    mom=result['momentum']['reference_common'];coverage=result['release_coverage']
    text+=f"""
## Momentum-transfer audit
For a constant spherical droplet, m=rho_l*pi*d^3/6 and A=pi*d^2/4. Newton's law gives m*dv/dt=(rho_g*C_D*A/2)*|w|*w + m*(1-rho_g/rho_l)*g, with w=u_g(r_d)-v_d. Dividing by mass gives 3*rho_g*C_D*|w|*w/(4*rho_l*d). Substituting C_D=24*F(Re_p)/Re_p and Re_p=rho_g*d*|w|/mu_g gives F(Re_p)*w/tau_p. At zero slip the implemented drag is exactly zero without evaluating C_D/|w|. The existing hybrid is retained: exact Stokes for Re_p<=0.01 and Schiller-Naumann above, with Re_p<800 checked at every RK stage. Documentation: https://doc.comsol.com/6.4/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.36.html . The small finite-Re correction below the declared Stokes threshold is neglected; this is explicitly a hybrid, not pure SN at every positive Re.

Current reference result is preserved, not fitted or tuned. Droplet mass {mom['droplet_mass_kg']:.12g} kg. Over the common interval, drag supplies axial impulse {mom['drag_impulse_y_p_N_s']:+.12g} N s, effective gravity supplies {mom['effective_gravity_impulse_y_p_N_s']:+.12g} N s, and droplet axial momentum changes {mom['momentum_change_y_p_kg_m_s']:+.12g} kg m/s. Equivalent velocity changes are {mom['drag_equivalent_delta_v_y_p_m_s']:+.9f} m/s from drag and {mom['gravity_equivalent_delta_v_y_p_m_s']:+.9f} m/s from effective gravity. RK4 impulse residual {mom['momentum_residual_RK4_y_p_N_s']:.6g} N s; independent saved-force trapezoid residual {mom['momentum_residual_trapezoid_y_p_N_s']:.6g} N s. The latter corresponds to {mom['residual_trapezoid_equivalent_delta_v_y_p_m_s']:.6g} m/s. Stage-quadrature closure is expected at roundoff because the same stages advance velocity; the independent trapezoid audit guards against treating that algebraic closure as independent validation. No residual is overwritten with zero.

The seven starts span {coverage['seed_span_mm']:.6f} mm. Measured entry-plane half-maximum crossings are x_p={coverage['left_halfmax_crossing_x_p_mm']:.6f} and {coverage['right_halfmax_crossing_x_p_mm']:.6f} mm (width {coverage['entry_FWHM_mm']:.6f} mm). The start span covers {100*coverage['seed_span_fraction_of_FWHM']:.2f}% of that width; no droplet-population weights are available. They are one reference and six unweighted sensitivity starts, not a measured whole-plume distribution. See initial_condition_coverage and release_coverage_profile.csv; physical asymmetry is retained.

The PZT-only plume can decelerate when aerodynamic slip opposes droplet motion and through gravity, but its surrounding gas was not measured to be stationary. Its spatial ensemble mean also reflects transport and sampling; it does not identify a unique individual drag history. Faster co-directed gas transfers momentum through drag alone. Finite response time leaves evolving slip; positive axial slip gives upward drag, but net upward acceleration requires that drag exceed effective weight. Smaller diameter reduces the Stokes response time quadratically and increases coupling; size is an assumed model parameter, not inferred by fitting this mean curve. Gravity contributes the signed impulse above. These are calculated conditional droplet paths in a prescribed mean gas field, not measured combined-flow results or automatically ensemble-mean droplet trajectories.

Analytical constant-gas Stokes position/velocity tests, zero-slip and sign tests, first-principles C_D force equivalence, instantaneous force balance, saved-force impulse quadrature, native support, inverse transforms and time-step refinement are recorded in validation/. Every trajectory CSV includes both drag/effective-gravity/total acceleration and force components, both momentum-budget components, response times, vector slip and particle Re. momentum_budget.csv contains full-path and common-interval budgets per droplet. No total plume momentum flux is inferred without droplet counts/mass loading.
"""
    (out/'README.md').write_text(text,encoding='utf-8')
    (out/'methods_results.tex').write_text(concise_latex(ctx,result),encoding='utf-8')
    def escape(text):return text.replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
    (out/'captions.tex').write_text('\n'.join('% '+stem+'\n\\caption{'+escape(caption)+'}\n' for stem,caption in captions.items()),encoding='utf-8')


def render_saved():
    """Regenerate presentation from saved, converged trajectories without reintegration."""
    from .aligned_codirected import load_context,prepare
    ctx=load_context();field,states,table,alignment,audit=prepare(ctx,write=False);out=ctx['out']
    result=json.loads((out/'results_summary.json').read_text());histories=[]
    for state in states:
        z=np.genfromtxt(out/'data'/f'release_{state["number"]:02d}.csv',names=True,delimiter=',')
        q={key:z[key] for key in z.dtype.names};nonpositive=np.flatnonzero(q['droplet_v_p_m_s']<=0);limit=int(nonpositive[0])+1 if len(nonpositive) else len(z);peak=int(np.argmax(q['y_p_mm'][:limit]))
        histories.append(dict(number=state['number'],reference=state['reference'],q=q,peak=peak))
    make_figures_and_report(ctx,field,histories,alignment,result)
    return result


def dimensionless_only():
    """Measured response/jet diagnostics without crossing or placement prerequisites."""
    from .aligned_codirected import load_context
    ctx=load_context();out=ctx['out'];(out/'data').mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    d=ctx['cfg']['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,ctx['props']));f=float(ctx['gas_cfg']['forcing_Hz']);D0=ctx['D0'];nu=ctx['props']['mu_air_Pa_s']/ctx['props']['rho_air_kg_m3'];valid=ctx['valid']&(ctx['uc']>0);uc=np.where(valid,ctx['uc'],np.nan)
    values=dict(ZNMF_downstream_mm=ctx['sx'],Uc_m_s=ctx['uc'],status=ctx['gas']['status'],supported_positive_centreline=valid,Re_D_c=uc*D0/nu,Sr_D_c=f*D0/uc,Stk_D_c=tau*uc/D0)
    write_csv(out/'data/jet_local_centreline_numbers_native.csv',values)
    write_csv(out/'data/response_scalars_independent_of_placement.csv',dict(nominal_estimated_diameter_um=[d*1e6],tau_stokes_s=[tau],Stk_f=[tau*f],forcing_frequency_Hz=[f],D0_m=[D0]))
    for key,label,stem in [('Re_D_c',r'$Re_{D,c}\ (-)$','jet_local_Reynolds'),('Sr_D_c',r'$Sr_{D,c}\ (-)$','jet_local_Strouhal')]:
        fig,ax=plt.subplots(figsize=(6.5,4.4),layout='constrained');ax.plot(ctx['sx'],values[key],color=RED);ax.set(xlabel=r'$X_j\ (\mathrm{mm})$',ylabel=label)
        for ext in ('pdf','svg','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig)
    return dict(tau_stokes_s=tau,Stk_f=tau*f,forcing_frequency_Hz=f,D0_m=D0,placement_required=False)


def concise_latex(ctx,result):
    """One self-contained scientific section, generated from calculated values."""
    entry=result['entry'];cross=result['crossing'];r=result['reference_comparison'];m=result['momentum']['reference_common']
    def sci(value):
        if value==0:return '0'
        exponent=int(np.floor(np.log10(abs(value))));return f'{value/10.**exponent:.7g}'+r'\times10^{'+str(exponent)+'}'
    text=r"""\section{Co-directed aerodynamic momentum transfer}
The separately measured PZT and ZNMF fields are paired as a prospective arrangement.
The original source pairing $X_j=y_{\rm lab}-\delta$ identifies the crossing at
$X_j=@CROSS@$ mm and the first supported native overtaking station at $X_j=@ENTRY@$ mm
($y_p=@YP@$ mm). One rigid placement aligns the jet centreline with the fitted plume axis
there. Both locally measured aligned droplet components initialise developed droplets;
neither the atomiser surface nor near-source velocities are relocated. Full-support
numerical fields retain their separate calibrations and masks; the report crop is unused.
Gravity is rotated with the saved frame. Interpolation requires four finite native corners
and trajectories stop at unsupported gas cells.

For a spherical, constant-diameter droplet,
$m_d=\rho_l\pi d^3/6$, $A_d=\pi d^2/4$ and
$\mathbf w=\mathbf u_g(\mathbf r_d)-\mathbf v_d$. Newton's law gives
\[
 \dot{\mathbf r}_d=\mathbf v_d,\qquad
 m_d\dot{\mathbf v}_d=\tfrac12\rho_g C_D A_d|\mathbf w|\mathbf w
 +m_d(1-\rho_g/\rho_l)\mathbf g_p.
\]
Using $Re_p=\rho_gd|\mathbf w|/\mu_g$ and $C_D=24F/Re_p$ yields the stable form
\[
 \dot{\mathbf v}_d=F(Re_p)\mathbf w/\tau_p+(1-\rho_g/\rho_l)\mathbf g_p,
 \qquad \tau_p=\rho_ld^2/(18\mu_g).
\]
The audited existing hybrid is retained: $F=1$ for $Re_p\le0.01$ and
$F=1+0.15Re_p^{0.687}$ above, with $Re_p<800$ checked at every RK stage.
The low-Re branch explicitly neglects the small finite-Re correction; zero slip has zero
drag without division by zero. No separate jet force, velocity addition or PZT-only
forcing is included. Nominal $d=@DIAM@\,\mu\mathrm m$ is an aperture/Rayleigh-based estimate,
not measured sizing. Saved water/air properties at $20\,{}^\circ\mathrm C$ and assumed
101.325 kPa are retained; diameter is constant and evaporation is neglected.

Over $y_p=@LOW@$--$@HIGH@$ mm, the reference velocity changes from @V0@ to @V1@ m/s
in @TIME@ ms, with transverse displacement @DX@ mm; these values are unchanged by the audit.
The single-droplet momentum budget is
\[
 m_d[\mathbf v_d(t)-\mathbf v_d(0)]=\int\mathbf F_D\,dt+
 \int\mathbf F_{g,\rm eff}\,dt+\boldsymbol\epsilon.
\]
Axial drag impulse is $@JD@\,\mathrm{N\,s}$, effective-gravity impulse
$@JG@\,\mathrm{N\,s}$ and momentum gain $@DP@\,\mathrm{kg\,m\,s^{-1}}$.
The RK4-stage residual is $@RK@\,\mathrm{N\,s}$; independent saved-force trapezoidal
quadrature gives $@TRAP@\,\mathrm{N\,s}$. Stage-quadrature closure is expected to
roundoff and is not independent physical validation. Zero-slip, force-direction,
analytical constant-gas Stokes position/velocity and step-refinement checks pass.

Faster gas transfers momentum through drag; finite response time leaves slip.
Upward drag produces net upward acceleration only when it exceeds effective weight.
Without the imposed jet, opposing slip and gravity can decelerate droplets, but the
PZT-only experiment establishes neither stationary surrounding air nor a unique drag
history from its spatial ensemble mean. Diameter strongly affects response through
$\tau_p\propto d^2$. The seven starts cover @COVER@\% of the measured entry-plane
half-maximum width and are unweighted sensitivity locations, not a whole-plume distribution.
$Stk_{D,\rm local}=\tau_p|\mathbf u_g|/D_0$ uses flow time $D_0/|\mathbf u_g|$;
$Stk_f=\tau_pf$. Local $Re_{D,c}=U_cD_0/\nu_g$ and $Sr_{D,c}=fD_0/U_c$ preserve gaps
and exclude nonpositive $U_c$, with documented $D_0=2$ mm and forcing $f=120$ Hz,
not the 10-Hz acquisition rate. These are not nozzle-exit quantities.
The prescribed gas is a steady, one-way-coupled measured mean: loading-induced gas changes
and phase-dependent motion are unresolved. Calculated paths are not measured combined-flow
results or automatically ensemble-mean droplet trajectories; no collection efficiency or
whole-plume momentum flux is established.
"""
    values={'DIAM':f"{ctx['cfg']['droplet_diameters_um'][0]:g}",'CROSS':f"{cross['ZNMF_downstream_mm']:.6f}",'ENTRY':f"{entry['ZNMF_downstream_mm']:.6f}",'YP':f"{entry['y_p_mm']:.6f}",'LOW':f"{r['lower_y_p_mm']:.6f}",'HIGH':f"{r['upper_y_p_mm']:.6f}",'V0':f"{r['entry_droplet_v_p_m_s']:.6f}",'V1':f"{r['end_droplet_v_p_m_s']:.6f}",'TIME':f"{r['transit_time_ms']:.6f}",'DX':f"{r['lateral_displacement_mm']:+.6f}",'JD':sci(m['drag_impulse_y_p_N_s']),'JG':sci(m['effective_gravity_impulse_y_p_N_s']),'DP':sci(m['momentum_change_y_p_kg_m_s']),'RK':sci(m['momentum_residual_RK4_y_p_N_s']),'TRAP':sci(m['momentum_residual_trapezoid_y_p_N_s']),'COVER':f"{100*result['release_coverage']['seed_span_fraction_of_FWHM']:.1f}"}
    for key,value in values.items():text=text.replace('@'+key+'@',value)
    return text
