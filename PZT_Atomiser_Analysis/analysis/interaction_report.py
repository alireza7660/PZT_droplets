"""Scientific report for interaction of developed droplets at the source-profile crossing."""
import json
import numpy as np


def write_report(out,cfg,result,stems):
    a=result['alignment'];audit=result['interaction_audit'];cross=audit['selected_crossing'];s=cross['ZNMF_downstream_mm'];y=cross['PZT_height_mm']
    initial=np.genfromtxt(out/'data/initial_conditions.csv',names=True,delimiter=',')
    centre=initial[initial['PZT_native_column']==a['PZT_reference_native_column']][0]
    intervals=np.genfromtxt(out/'data/common_height_comparison.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')
    reference_index=list(result['release_columns']).index(a['PZT_reference_native_column'])
    reference=intervals[reference_index]
    lo,hi=result['common_interval_y_mm'];tlo,thi=result['common_travel_time_range_ms'];dvlo,dvhi=result['common_delta_v_range_m_s']
    readme=f'''# Developed-droplet interaction with the upward measured ZNMF jet

## Source tracing and the correction

The archived comparison script and its plotted CSVs are traced in
interaction_plane_audit.json. Its curves use ZNMF downstream distance s and
PZT plume height y as separate source coordinates on an overlaid 0..27 mm axis.
The reproduced curve samples match those archived CSVs exactly (75 gas and
315 PZT points); PZT is red in the current presentation. Velocities are already
SI and are never recalibrated. The reproduction is source_comparison_reproduced.

Adjacent VALID source intervals give one gas-overtaking crossing:
ZNMF s={s:.12g} mm, PZT y={y:.12g} mm, axial velocity
{cross['common_velocity_m_s']:.12g} m/s. No 10 mm constant or near-atomiser
velocity is inserted. Roots are solved on the union of native profile knots,
with no interpolation across invalid centreline or PZT gaps. The first supported
gas-overtaking root is selected; all roots and brackets are exported.

The superseded prediction instead paired the lowest measurement planes and
sampled a near-entry gas reversal region. That was a different initial condition
and registration, not the requested interaction of already-developed droplets.
Its generated predictions and code/configuration snapshots are archived under
{cfg['co_directed']['archive_manifest']}. Original source experiments remain unchanged.

## Frames and alignment assumption

The common frame keeps PZT x rightward and y upward. Gas source Xj is downstream,
Yj is transverse relative to its documented geometric nozzle centre. At interaction,
the measured gas centreline Yj,c={a['ZNMF_interaction_centreline_Yj_mm']:.9g} mm
is aligned to the documented fixed PZT reference x={a['plume_reference_x_mm']:.8g} mm.
Positions: x=x_ref+Yj,c-Yj; y=y_PZT,int+Xj-s_ZNMF,int.
Vectors: u_g=-V_ZNMF, v_g=+U_ZNMF (+90 degree rotation).
The numerical pairing of the independent axial coordinates is an explicit
comparison-based registration assumption, not evidence of coincident physical
exits. The configurable axial offset is {a['relative_axial_offset_mm']:.8g} mm.
Centreline alignment is applied once at the interaction plane; no downstream
straightening or ongoing particle recentering is applied. Native transverse offsets
are retained. All transforms act on analysis copies; source values/masks do not change.

## Initial conditions exported BEFORE integration

data/initial_conditions.csv gives source coordinates, PZT rows/columns and height
weights, gas-cell indices/weights, common-frame position, measured droplet components,
sampled gas components, vector slip, Re and initial drag/net acceleration.
data/initial_interpolation_stencils.json contains the original corner values.
The exact crossing falls between PZT rows {int(centre['PZT_lower_native_row'])}
and {int(centre['PZT_upper_native_row'])}; both components are linearly interpolated
in height at each native column, with no adjacent-column averaging. These are
developed droplets at the interaction height, not droplets reset to their initial
near-atomiser velocity or relocated atomiser exits. Vector sample counts are not
interpreted as droplet concentration.

Reference release: droplet (u,v)=({centre['measured_droplet_u_m_s']:.9g},
{centre['measured_droplet_v_m_s']:.9g}) m/s, sampled gas (u,v)=
({centre['sampled_gas_u_m_s']:.9g},{centre['sampled_gas_v_m_s']:.9g}) m/s.
The plotted gas centreline is {centre['plotted_centreline_Uc_m_s']:.9g} m/s;
sampled-minus-plotted axial difference is {centre['sampled_minus_plotted_axial_m_s']:.9g} m/s.
Initial net axial acceleration is {centre['initial_net_axial_acceleration_m_s2']:.9g} m/s2.
The crossing equates axial profiles, not vector velocities. In ZNMF pipeline.py,
extract_centreline finds subgrid peak positions and linearly interpolates velocities
there. Interpolating that moving-centreline curve downstream differs from bilinear
sampling at a single interpolated physical transverse position. Off-axis releases
add real transverse-profile differences. No gas value is changed to force equality.

## Actual predictions and common supported comparison

Completed cases: {result['successfully_integrated_cases']}. Terminations:
{result['termination_reasons']}. All histories begin at the interaction plane.
The common first-ascent interval is y={lo:.9g}..{hi:.9g} mm, limited by BOTH
measured PZT profile support and all computed ascending trajectories. Velocity
changes are {dvlo:.9g}..{dvhi:.9g} m/s and travel times {tlo:.9g}..{thi:.9g} ms.
Reference release: delta_v={reference['delta_v_m_s']:.9g} m/s over
{reference['travel_time_s']*1e3:.9g} ms, ending v={reference['final_v_m_s']:.9g} m/s.
These compare the same height bounds, not different endpoint displacements.
The full case data and common-height table retain individual results and termination
times. Interpolation to common diagnostic bounds stays within computed support.
The measured PZT-only red curve is a fixed-line ensemble mean, not a tracked
particle trajectory; its height trend can also reflect lateral plume displacement.

## Model and limits (prescribed gas approximation)

The measured two-component ZNMF mean alone prescribes gas velocity. It neglects
plume-induced gas changes, jet unsteadiness, evaporation and droplet-size evolution.
PZT velocities initialise particles once; velocities are not added and no downstream
PZT-only profile is imposed. Nominal d=10 um is a model-based Rayleigh/aperture
estimate (5 um aperture, assumed equal jet diameter, 1.89 factor ->9.45 um rounded),
not measured sizing. Existing 20 C water/air properties are retained; pressure
101325 Pa and equal water temperature are assumptions, with dry-air properties
approximating unknown humidity.

The SI force model uses vector slip, gravity and buoyancy with spherical drag:
a=-g_eff*e_y+slip*f(Re)/tau_Stokes. Re=rho_air*d*norm(slip)/mu_air;
f=1 at Re<=0.01, otherwise 1+0.15*Re^0.687, restricted to Re<800. Gas interpolation
requires four finite paired native corner vectors; no gaps/boundaries are extrapolated.
Predictions stop at unsupported gas support or the configured time horizon.
Adaptive numerical time-step refinement enforces the existing convergence tolerances
without changing the force law or initial measurements.

tau_Stokes=rho_water*d^2/(18*mu_air) is a Stokes-limit baseline. The exported
tau_Stokes/f(Re) is a drag/slip coefficient timescale, not an exact finite-Re relaxation
time. Stk_f=tau_Stokes*f uses documented 120 Hz jet actuation (no 2*pi), not
10 Hz acquisition or 190 kHz atomiser frequency. Local trajectory Stk_D uses
D0/norm(local gas vector) with D0=2 mm. Local gas Re_D,c=Uc*D0/nu and
Sr_D,c=f*D0/Uc use positive finite measured Uc with gaps retained; they are not
nozzle-exit Re0/Sr0. Stokes Stk and Strouhal Sr are distinct. Scalars are tabulated.
Supporting PZT-only particle Re has an explicit zero-background-gas reference,
not a measurement of ambient PZT air. No collection efficiency, lung deposition,
or experimental validation of the combined flow is inferred. The remaining axial
and transverse registration assumptions and nominal diameter limit physical certainty.

## Reproduce and files

From the PZT_droplets workspace root:
.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook

report_manifest.json explicitly selects active PDF/600-dpi PNG figures and data.
source_manifest.json records original paths/configuration, transformations and hashes.
interaction_plane_audit.json records source reproduction, the computed root, and
pre-integration checks. data/ contains initial audits, full trajectories, the exact
main-comparison sources, common-height comparisons, force intervals, scalars and
local dimensionless diagnostics. validation/ contains convergence/preservation checks.
Archived H/L/C and superseded near-atomiser predictions are excluded from the notebook.
'''
    (out/'README.md').write_text(readme,encoding='utf-8')
    tex=rf'''\section{{Interaction of developed PZT droplets with the measured upward jet}}
The archived comparison is reproduced from its authoritative PZT fixed-line and
ZNMF measured-centreline sources, using their original independent axial coordinates.
Supported piecewise-linear interpolation gives a unique gas-overtaking crossing at
$s_{{ZNMF}}={s:.9g}\,$mm and $y_{{PZT}}={y:.9g}\,$mm, with common axial velocity
${cross['common_velocity_m_s']:.9g}\,\mathrm{{m\,s^{{-1}}}}$. No approximate 10 mm value is imposed.
Gas downstream and PZT height coordinates are paired as in that comparison; this
is an explicit registration assumption, not evidence of coincident physical exits.
The measured gas centreline at the interaction plane is aligned once to the
documented PZT fixed-line reference. Common positions and vectors obey
\begin{{align}}
x&=x_{{PZT,ref}}+Y_{{j,c,int}}-Y_j,\qquad
y=y_{{PZT,int}}+X_j-s_{{ZNMF,int}},\\
u_g&=-V_{{ZNMF}},\qquad v_g=U_{{ZNMF}}.
\end{{align}}
Transverse offsets and downstream centreline wandering remain unchanged. Saved
SI velocities, signs, masks and sample populations are preserved without recalibration.

Seven representative initial states use both measured PZT components at the
interaction height, interpolating only between adjacent supported native height
planes at each selected column. No near-atomiser velocity is reused. Before
integration, source indices/weights, transformed positions, gas samples, vector
slip and acceleration are exported. For the reference release, droplet axial
velocity is {centre['measured_droplet_v_m_s']:.8g} m/s and sampled gas axial velocity is
{centre['sampled_gas_v_m_s']:.8g} m/s, compared with plotted centreline
{centre['plotted_centreline_Uc_m_s']:.8g} m/s. The two interpolation operations differ
when the measured centreline moves between native columns. Transverse velocities
also differ: droplet {centre['measured_droplet_u_m_s']:.8g} m/s and gas
{centre['sampled_gas_u_m_s']:.8g} m/s. Thus an axial-profile crossing is not zero vector slip.
The reference release initially has net axial acceleration
{centre['initial_net_axial_acceleration_m_s2']:.8g} m/s$^2$; acceleration is not forced positive.

Nominal $d=10\,\mu$m is a configurable model estimate from an assumed aperture-sized
jet ($5\,\mu$m) and Rayleigh factor 1.89, rounded from $9.45\,\mu$m, not measured
sizing. Existing 20-degree-C properties are used. The prescribed two-component
measured mean gas neglects plume-induced gas changes, unsteadiness, evaporation
and size evolution. PZT velocities initialise the particles once and are never
added to the gas or imposed downstream. The equations are
\begin{{align}}
\dot{{\boldsymbol q}}&=\boldsymbol v_d,\\
\dot{{\boldsymbol v}}_d&=-g(1-\rho_a/\rho_w)\boldsymbol e_y+
\frac{{f(Re_p)}}{{\tau_p}}(\boldsymbol u_g-\boldsymbol v_d),\\
Re_p&=\rho_a d\lVert\boldsymbol u_g-\boldsymbol v_d\rVert/\mu_a,
\qquad\tau_p=\rho_wd^2/(18\mu_a).
\end{{align}}
Use $f=1$ for $Re_p\leq0.01$ and $f=1+0.15Re_p^{{0.687}}$ otherwise,
with $Re_p<800$. Only four-finite-vector gas cells are interpolated; boundaries
and invalid regions are never extrapolated. Gravity/buoyancy remain downward
in the upward-positive common frame. Numerical refinement is checked independently
without tuning physical inputs.

All {result['successfully_integrated_cases']} integrations complete. Over the SAME supported first-ascent
interval $y={lo:.8g}$--${hi:.8g}\,$mm, axial velocity changes span
{dvlo:.8g}--{dvhi:.8g} m/s and travel times span {tlo:.8g}--{thi:.8g} ms.
For the reference release, the change is {reference['delta_v_m_s']:.8g} m/s over
{reference['travel_time_s']*1e3:.8g} ms. These are common-height comparisons, not
performance rankings of displacements at different boundaries/times. The measured
red PZT-only fixed-line ensemble profile remains an experimental baseline, not
a tracked individual trajectory or tracked centreline. Differences between it
and predictions do not validate the combined flow or identify a measured size.

The scalar Stokes-limit time and $Stk_f=\tau_p f$ are tabulated, using documented
ZNMF actuation $f=120\,$Hz rather than acquisition or atomiser frequency.
Local trajectory $Stk_D=\tau_p\lVert\boldsymbol u_g\rVert/D_0$ uses the local
nozzle-length gas timescale. $\tau_p/f(Re_p)$ is a coefficient timescale, not an
exact finite-Re relaxation time. Local gas $Re_{{D,c}}=U_cD_0/\nu_a$ and
$Sr_{{D,c}}=fD_0/U_c$ use positive finite measured $U_c$, with unsupported gaps
retained and $D_0=2\,$mm; they are not nozzle-exit $Re_0/Sr_0$.
Stk and Sr are distinct. No collection efficiency, lung deposition or
experimental validation is claimed. The comparison-coordinate registration,
single-plane centreline alignment, prescribed gas and estimated diameter remain limitations.
'''
    (out/'methods_results.tex').write_text(tex,encoding='utf-8')
    captions={
        'source_comparison_reproduced':'Measured source profiles reproduced from the actual archived-comparison data selections, with PZT in red. The independent gas downstream and PZT plume coordinates are overlaid for this comparison; the crossing is computed only within adjacent valid intervals, not imposed at 10 mm.',
        'measured_and_predicted_axial_comparison':'Measured ZNMF centreline velocity (black dashed), measured PZT-only fixed-line ensemble velocity (red), predicted droplet velocity for the reference release (blue), and gas sampled on that same trajectory (purple dotted). Pale blue curves show other representative releases. The vertical line marks interaction. Measured gas centreline values and sampled gas are distinct. Predictions begin with developed PZT velocities at that plane; no upstream PZT profile is imposed on them.',
        'trajectories_upward_jet':'Conditional developed-droplet trajectories over signed upward-oriented measured gas axial contours. The measured centrelines are aligned once at the numerical interaction plane; transverse offsets and subsequent wandering remain. No experimental gas or droplet velocity is tuned or added.',
        'trajectories_entry_detail':'Detail of the first 3 mm downstream of interaction. Initial positions are identified by physical x coordinate; measured developed PZT velocity pairs initialise these releases. The plot is not a relocated atomiser launch.',
        'axial_velocity_along_trajectories':'Axial droplet and gas velocities sampled at the same evolving trajectory positions. Particle ascent/descent use solid/dashed branches; gas traces use dotted/dash-dot branches in matching colours. They are distinct from measured centreline peak values.',
        'axial_slip_and_acceleration':'Axial vector-slip component and net axial acceleration (solid/ascent, dashed/descent), with drag alone dotted. The positive black threshold is drag equal to buoyancy-corrected gravity, not upward gravity. Upward drag need not imply net upward acceleration. Acceleration uses a symmetric-log ordinate, linear within plus/minus 20 m/s squared, retaining all signs.',
        'particle_Re_and_local_Stokes':'Slip-based particle Reynolds number and supporting local Stokes diagnostic along integrated trajectories. The latter uses the Stokes-limit time and local gas-vector magnitude over the 2 mm nozzle length; it is not capture or transport efficiency.',
        'PZT_only_measured_axial_profile':'Measured PZT-only fixed-line ensemble-mean axial velocity in red. This is an experimental spatial baseline, not an individual trajectory; lateral plume displacement can affect the fixed-line trend.',
        'jet_local_dimensionless_numbers':'Local measured-centreline gas Reynolds and Strouhal numbers and supporting particle Stokes diagnostic. The axial coordinate is transformed using the documented interaction alignment. Positive finite measured Uc is used with gaps retained; no nozzle-exit velocity or Re0/Sr0 is inferred. Scalars are tabulated, not constant-line plots.'}
    (out/'captions.tex').write_text(''.join('% '+stem+'\n'+r'\caption{'+captions[stem]+'}\n\n' for stem in stems),encoding='utf-8')
