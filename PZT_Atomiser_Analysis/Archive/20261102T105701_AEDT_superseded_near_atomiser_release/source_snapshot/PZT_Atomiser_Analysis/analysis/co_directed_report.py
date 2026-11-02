"""Dimensionless tables and concise scientific text for the co-directed package."""
import json
import numpy as np
from .combined_flow_prediction import ROOT,write_csv,write_json,sha
from .droplet_dynamics import relaxation_time,reynolds
from .plotting import plt


def dimensionless(out,cfg,p,D0,alignment):
    config_path=(ROOT/cfg['znmf_configuration']).resolve();gcfg=json.loads(config_path.read_text())
    source=config_path.parent.parent/'tables/FOV_1_2_3_centreline_and_width.csv'
    z=np.genfromtxt(source,names=True,delimiter=',',dtype=None,encoding='utf-8')
    valid=(z['status']=='measured')&np.isfinite(z['Uc_m_per_s'])&(z['Uc_m_per_s']>0)
    uc=np.where(valid,z['Uc_m_per_s'],np.nan);f=float(gcfg['forcing_Hz'])
    d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,p));nu=p['mu_air_Pa_s']/p['rho_air_kg_m3']
    y=z['x_over_D0']*D0*1e3+alignment['translation_m'][1]*1e3
    values=dict(x_ZNMF_over_D0=z['x_over_D0'],aligned_y_mm=y,Uc_measured_m_s=z['Uc_m_per_s'],supported=valid,
                Re_D_c=uc*D0/nu,Sr_D_c=f*D0/uc,Stk_D_c=tau*uc/D0,status=z['status'])
    write_csv(out/'data/jet_local_centreline_numbers.csv',values)
    scalars=dict(nominal_diameter_um=[d*1e6],tau_Stokes_s=[tau],tau_Stokes_ms=[tau*1e3],
        jet_actuation_frequency_Hz=[f],acquisition_frequency_Hz_not_used=[gcfg['acquisition_Hz']],D0_m=[D0],
        Stk_f=[tau*f],rho_air_kg_m3=[p['rho_air_kg_m3']],rho_water_kg_m3=[p['rho_water_kg_m3']],
        mu_air_Pa_s=[p['mu_air_Pa_s']],nu_air_m2_s=[nu])
    write_csv(out/'data/scalar_quantities.csv',scalars)
    # Retain the PZT-only diagnostic as explicit supporting reference data, not gas measurement.
    with np.load(ROOT/cfg['pzt_statistics']) as a:
        re=reynolds(d,np.stack([a['mean_u'],a['mean_v']],axis=-1),[0,0],p)
        write_csv(out/'data/supporting_PZT_zero_background_Re.csv',dict(x_mm=a['x_mm'].ravel(),y_mm=a['y_mm'].ravel(),
            mean_u_m_s=a['mean_u'].ravel(),mean_v_m_s=a['mean_v'].ravel(),Re_p_zero_background_reference=re.ravel()))
    fig,axes=plt.subplots(3,1,figsize=(6.6,7.4),sharex=True,layout='constrained')
    for ax,key,label in zip(axes,['Re_D_c','Sr_D_c','Stk_D_c'],[r'$Re_{D,c}\ (-)$',r'$Sr_{D,c}\ (-)$',r'$Stk_{D,c}\ (-)$']):
        ax.plot(y,values[key],color='tab:blue');ax.set_ylabel(label)
    axes[-1].set_xlabel(r'$y\ (\mathrm{mm})$');stem='jet_local_dimensionless_numbers'
    for ext in ('pdf','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)
    return dict(tau_Stokes_s=tau,Stk_f=tau*f,actuation_frequency_Hz=f,D0_m=D0,nu_air_m2_s=nu,
        positive_measured_centreline_rows=int(valid.sum()),unsupported_rows=int((~valid).sum()),
        source=str(source),source_sha256=sha(source),definitions='Re_D,c=Uc D0/nu; Sr_D,c=f D0/Uc; Stk_D,c=tau_Stokes Uc/D0; trajectory Stk_D=tau_Stokes norm(gas vector)/D0'),stem


def write_report(out,cfg,result,stems):
    a=result['alignment'];n=result['dimensionless'];interval=result['common_interval_y_mm'];dv=result['common_delta_v_range_m_s'];travel=result['common_travel_time_range_ms']
    readme=f'''# Upward co-directed prescribed-field prediction

Actual completed integrations: {result['successfully_integrated_cases']}. The experiments
were separate; no combined-flow validation, collection efficiency or lung-deposition
claim is made. The obsolete horizontal side-jet study is not regenerated.

## Geometry and measured initial states

Prediction x is horizontal/rightward and y upward. The documented PZT fixed-line
reference is x={a['plume_reference_x_mm']:.8g} mm, native column {a['PZT_reference_native_column']},
through the global mean-v maximum (methods/vertical_profile_selection.json in the
authoritative PZT results). It is an empirical plume reference, not a located outlet.
Nominal gas axis is aligned to it; real gas-centreline wandering is not straightened.
Gas rotation is +90 degrees: (x,y)=(x_ref-Yj, Xj+translation_y),
(u_g,v_g)=(-V_ZNMF,+U_ZNMF). No source velocity is scaled or sign-tuned.
The PZT lowest measurement plane y={a['PZT_lowest_plane_y_mm']:.8g} mm and first
supported gas downstream plane Xj={a['ZNMF_first_supported_Xj_mm']:.8g} mm are aligned
with axial offset {a['relative_axial_offset_mm']:.8g} mm. This is an assumed coincidence
of MEASUREMENT PLANES, not physical exits. Offset stays configurable in co_directed.
All transformations apply to analysis copies. Exact cardinal rotations avoid
floating-point boundary artefacts; no measured-field extension is introduced.

Seven native PZT columns {result['release_columns']} at row {result['release_row']} supply
the first common supported release plane, with both signed measured components.
All starts and source indices are in data/trajectory_summary.csv. These are
representative measured velocity samples, not tracked droplet identities or
droplet concentration. Vector sample counts are never used as number concentration.
No unmeasured travel from the atomiser surface is simulated. PZT velocities are
used once for initialisation and are never added to gas velocities or imposed again.

## Actual results

Initial sampled gas axial velocities range {result['initial_gas_axial_velocity_range_m_s'][0]:.7g}
to {result['initial_gas_axial_velocity_range_m_s'][1]:.7g} m/s. Negative local axial gas
measurements near the nominal axis at the first plane are preserved, despite the
nominal upward jet direction. They are not masked or made positive to force coflow.
Cases ending with downward droplet velocity: {result['reversed_cases']}.
Cases with any upward drag: {result['cases_with_upward_drag']}; cases with any net
upward axial acceleration: {result['cases_with_net_upward_acceleration']}.
Maximum heights span {result['peak_height_range_mm'][0]:.8g} to
{result['peak_height_range_mm'][1]:.8g} mm. Terminations: {result['termination_reasons']}.

Comparisons use the SAME first-ascent height interval y={interval[0]:.8g} to
{interval[1]:.8g} mm, not differing endpoints. Travel times are {travel[0]:.7g} to
{travel[1]:.7g} ms; axial velocity changes are {dv[0]:.7g} to {dv[1]:.7g} m/s.
See data/common_height_comparison.csv for each case. Values at common bounds use
interpolation along computed first-ascent histories, never gas-field extrapolation.
Full curves retain later descent as separate dashed branches; they are not treated
as a single-valued tracked velocity-height profile. The red experimental PZT-only
ensemble profile is separate and is not interpreted as an individual trajectory.

Drag acts on vector slip. Positive axial slip means upward drag, but net upward
acceleration requires a_D,y > g*(1-rho_air/rho_water). The CSVs distinguish drag
and net acceleration, including upward-drag/net-downward intervals. Positive net
acceleration while a particle is descending means its downward motion decelerates;
it does not establish upward transport. Supporting force intervals are tabulated.

## Model, properties and dimensionless definitions

This prescribed mean-gas-field model neglects plume-induced gas changes, jet
unsteadiness, evaporation and droplet-size evolution. Spherical rigid-interface
drag and one-way coupling are assumptions. Diameter is nominal 10 um, motivated
by aperture 5 um times Rayleigh factor 1.89 = 9.45 um, rounded to 10 um; it is a
model estimate, not measured sizing or a validated distribution. Water and air
properties are the existing saved 20 C values. Pressure 101325 Pa and equal water
temperature are assumptions; dry-air properties approximate unknown humidity.

Acceleration is -g_eff e_y + slip*f(Re_p)/tau_Stokes. Use f=1 for Re<=0.01 and
f=1+0.15 Re^0.687 otherwise, with Re<800 enforced. Zero slip is finite. The
gas interpolant uses only cells with four finite paired source vectors and never
extrapolates through gaps or boundaries. Integration continues through turns but
stops at unsupported gas support, the drag validity limit or configured time horizon.

Re_p=rho_air*d*norm(vector slip)/mu_air. tau_Stokes=rho_water*d^2/(18*mu_air)
is the Stokes-limit response time; tau_drag_coefficient=tau_Stokes/f(Re) is only
an instantaneous drag/slip coefficient timescale, not an exact finite-Re relaxation
time. Local trajectory Stk_D=tau_Stokes*norm(gas velocity)/D0 uses the chosen
timescale D0/norm(gas velocity). At zero gas magnitude this timescale is infinite
and Stk_D is zero. Stk_f=tau_Stokes*f={n['Stk_f']:.9g} uses documented jet actuation
f={n['actuation_frequency_Hz']:g} Hz, not acquisition 10 Hz or atomiser 190 kHz, and
does not include 2*pi. Scalars are in data/scalar_quantities.csv, not constant-line plots.
Local gas Re_D,c=Uc*D0/nu and Sr_D,c=f*D0/Uc use D0=2 mm and positive finite
measured Uc, retaining gaps. Stk denotes Stokes number; Sr denotes Strouhal number.
No nozzle-exit Re0 or Sr0 is inferred. Supporting PZT-only Re uses an explicit
ZERO-BACKGROUND-GAS REFERENCE, not a measurement of ambient PZT gas velocity.

## Reproduce and report

From PZT_droplets workspace root:
`.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook`

Active report selection is report_manifest.json. figures/ contains seven title-free
PDF and 600-dpi PNG sets. methods_results.tex and captions.tex are report snippets.
data/ contains seven histories, trajectory summary, common-interval comparison,
force intervals, scalar table, local gas numbers, PZT-only baseline and explicitly
conditional PZT reference Re. Validation records cover common-frame rotation,
force signs, support, step convergence and unchanged source hashes. The archive
manifest is {cfg['co_directed']['archive_manifest']}; old H/L/C identifiers and
constant-line figures are excluded. Source paths, processing configurations,
original frame definitions and applied analysis-copy transforms are in source_manifest.json.
'''
    (out/'README.md').write_text(readme,encoding='utf-8')
    tex=rf'''\section{{Co-directed response of nominal PZT droplets to the measured ZNMF jet}}
The experiments were separate. Final SI PZT droplet and ZNMF gas velocities are
used without recalibration, rescaling, sign tuning or addition of the fields.
The ZNMF downstream direction is rotated upward: $u_g=-V$, $v_g=U$;
$x=x_{{ref}}-Y_j$, $y=X_j+y_{{PZT,low}}-X_{{j,first}}+\delta_y$.
The nominal axis is aligned to the documented empirical PZT reference
$x_{{ref}}={a['plume_reference_x_mm']:.8g}\,$mm without straightening measured wandering.
Because the PZT physical outlet is not established, the lowest PZT measurement
plane ({a['PZT_lowest_plane_y_mm']:.8g} mm) and first supported gas plane
({a['ZNMF_first_supported_Xj_mm']:.8g} mm in the gas frame) are made coincident,
with configurable $\delta_y={a['relative_axial_offset_mm']:.6g}\,$mm. This does not align physical exits.
Seven supported native PZT positions across the plume initialise both measured
velocity components once; these are not individually tracked droplets or number
concentration measurements. No travel from the atomiser surface is inferred.

For nominal $d=10\,\mu$m, motivated by $5\,\mu$m aperture times 1.89 and
rounding the $9.45\,\mu$m estimate, the equations are
\begin{{align}}
\dot{{\boldsymbol q}}&=\boldsymbol v_d,\\
\dot{{\boldsymbol v}}_d&=-g(1-\rho_a/\rho_w)\boldsymbol e_y+
\frac{{f(Re_p)}}{{\tau_p}}(\boldsymbol u_g-\boldsymbol v_d),\\
Re_p&=\rho_a d\lVert\boldsymbol u_g-\boldsymbol v_d\rVert/\mu_a,
\qquad\tau_p=\rho_w d^2/(18\mu_a).
\end{{align}}
The diameter is a configurable model estimate, not measured sizing. Existing
20-degree-C properties are retained. The drag multiplier is unity for $Re_p\leq0.01$
and $1+0.15Re_p^{{0.687}}$ otherwise, with $Re_p<800$. Gravity and buoyancy
remain in the upward-positive physical frame. Gas interpolation requires four
finite paired vectors per source cell; boundaries and gaps are never extrapolated.
This prescribed-mean-gas approximation neglects plume-induced gas changes, jet
unsteadiness, evaporation and droplet-size evolution. The combined flow is not
experimentally validated.

All {result['successfully_integrated_cases']} integrations complete. Initial sampled gas axial velocities
span {result['initial_gas_axial_velocity_range_m_s'][0]:.5g} to {result['initial_gas_axial_velocity_range_m_s'][1]:.5g} m/s;
the saved near-entry signed reversal region is retained, despite the nominal upward axis.
{result['reversed_cases']} trajectories end with downward droplet velocity. Maximum heights
span {result['peak_height_range_mm'][0]:.6g}--{result['peak_height_range_mm'][1]:.6g} mm.
Cases exhibiting upward drag: {result['cases_with_upward_drag']}; cases exhibiting
net upward axial acceleration: {result['cases_with_net_upward_acceleration']}.
Over the common first-ascent interval $y={interval[0]:.6g}$--${interval[1]:.6g}\,$mm,
travel times span {travel[0]:.6g}--{travel[1]:.6g} ms and axial velocity changes
span {dv[0]:.6g}--{dv[1]:.6g} m/s. These are common-interval comparisons, not
performance comparisons of differently terminated endpoints. The full histories
retain subsequent descent; branch styles distinguish ascent and descent.

Upward drag requires positive axial slip; net upward acceleration additionally
requires $a_{{D,y}}>g(1-\rho_a/\rho_w)$. Force-sign intervals are tabulated.
Net upward acceleration of a descending particle can mean deceleration of its
downward motion, not upward transport. The measured fixed-line PZT-only mean
profile is shown separately in red and is not a tracked single-droplet trajectory.

$\tau_p={n['tau_Stokes_s']*1e3:.6g}\,$ms is a Stokes-limit time; the exported
$\tau_p/f(Re_p)$ is a finite-Re drag/slip coefficient timescale, not an exact
relaxation time. Scalars are tabulated. $Stk_f=\tau_p f={n['Stk_f']:.7g}$ uses
ZNMF actuation $f={n['actuation_frequency_Hz']:g}\,$Hz without $2\pi$; neither acquisition
nor atomiser frequency is substituted. Local trajectory $Stk_D=\tau_p\lVert\boldsymbol u_g\rVert/D_0$
uses $D_0/\lVert\boldsymbol u_g\rVert$. Gas local $Re_{{D,c}}=U_cD_0/\nu_a$
and $Sr_{{D,c}}=fD_0/U_c$ use positive finite measured $U_c$ with retained gaps,
$D_0=2\,$mm, and do not imply nozzle-exit $Re_0$ or $Sr_0$. Stk and Sr are
distinct. Supporting PZT-only Re is conditional on explicitly assumed zero
background gas, not a measured ambient-gas condition. No collection efficiency,
lung deposition or combined-experiment validation is claimed.
'''
    (out/'methods_results.tex').write_text(tex,encoding='utf-8')
    captions={
      'trajectories_upward_jet':'Predicted trajectories over the signed upward-oriented measured gas axial component. Nominal axes and measurement reference planes are aligned; neither local negative gas values nor measured centreline wandering are removed. The complete measured gas support is shown; near-entry paths are resolved separately.',
      'trajectories_entry_detail':'Entry-region detail of the same trajectories and gas contours. Markers show measured PZT initial positions. These are conditional representative-droplet predictions; negative near-entry gas velocities and trajectory reversals are retained. Initial horizontal positions are identified in millimetres.',
      'axial_velocity_along_trajectories':'Droplet axial velocity and gas axial velocity sampled on the same trajectories versus height. Particle ascent/descent are solid/dashed; corresponding gas traces are dotted/dash-dot in the same colour. This is a time-parametrised multi-branch comparison, not a single-valued experimental trajectory profile.',
      'axial_slip_and_acceleration':'Signed axial slip and acceleration along predicted histories. Net axial acceleration uses solid/ascent and dashed/descent branches; faint dotted traces show drag alone. The black dotted positive reference is the drag threshold a_D,y=g_eff for zero net acceleration, not an upward gravitational force. Upward drag is not sufficient for net upward acceleration. Acceleration uses a symmetric-log ordinate with a linear region within plus/minus 20 m/s squared; signed values are not clipped.',
      'particle_Re_and_local_Stokes':'Evolving in-plane-slip particle Reynolds number and supporting local Stokes number. The latter uses the Stokes-limit time and local gas-vector magnitude with D0=2 mm; it is not a transport or capture criterion. Branch styles distinguish ascent and descent.',
      'PZT_only_measured_axial_profile':'Measured PZT-only ensemble mean axial velocity along the existing fixed vertical line, shown in red. The jet was absent and ambient gas velocity was not measured. This is not a tracked droplet trajectory or tracked plume centreline.',
      'jet_local_dimensionless_numbers':'Local measured-centreline jet Reynolds and Strouhal numbers and supporting particle Stokes comparison versus aligned height. Only positive finite measured centreline values are used, with unsupported gaps retained. No nozzle-exit velocity or Re0/Sr0 is inferred. Constant scalar quantities are tabulated separately.'}
    (out/'captions.tex').write_text(''.join('% '+s+'\n'+r'\caption{'+captions[s]+'}\n\n' for s in stems),encoding='utf-8')
