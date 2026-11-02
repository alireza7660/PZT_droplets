# Developed-droplet interaction with the upward measured ZNMF jet

## Source tracing and the correction

The archived comparison script and its plotted CSVs are traced in
interaction_plane_audit.json. Its curves use ZNMF downstream distance s and
PZT plume height y as separate source coordinates on an overlaid 0..27 mm axis.
The reproduced curve samples match those archived CSVs exactly (75 gas and
315 PZT points); PZT is red in the current presentation. Velocities are already
SI and are never recalibrated. The reproduction is source_comparison_reproduced.

Adjacent VALID source intervals give one gas-overtaking crossing:
ZNMF s=10.0860743924 mm, PZT y=10.0860743924 mm, axial velocity
0.964010941135 m/s. No 10 mm constant or near-atomiser
velocity is inserted. Roots are solved on the union of native profile knots,
with no interpolation across invalid centreline or PZT gaps. The first supported
gas-overtaking root is selected; all roots and brackets are exported.

The superseded prediction instead paired the lowest measurement planes and
sampled a near-entry gas reversal region. That was a different initial condition
and registration, not the requested interaction of already-developed droplets.
Its generated predictions and code/configuration snapshots are archived under
PZT_Atomiser_Analysis\Archive\20261102T105701_AEDT_superseded_near_atomiser_release\archive_manifest.json. Original source experiments remain unchanged.

## Frames and alignment assumption

The common frame keeps PZT x rightward and y upward. Gas source Xj is downstream,
Yj is transverse relative to its documented geometric nozzle centre. At interaction,
the measured gas centreline Yj,c=-0.258909028 mm
is aligned to the documented fixed PZT reference x=14.54665 mm.
Positions: x=x_ref+Yj,c-Yj; y=y_PZT,int+Xj-s_ZNMF,int.
Vectors: u_g=-V_ZNMF, v_g=+U_ZNMF (+90 degree rotation).
The numerical pairing of the independent axial coordinates is an explicit
comparison-based registration assumption, not evidence of coincident physical
exits. The configurable axial offset is 0 mm.
Centreline alignment is applied once at the interaction plane; no downstream
straightening or ongoing particle recentering is applied. Native transverse offsets
are retained. All transforms act on analysis copies; source values/masks do not change.

## Initial conditions exported BEFORE integration

data/initial_conditions.csv gives source coordinates, PZT rows/columns and height
weights, gas-cell indices/weights, common-frame position, measured droplet components,
sampled gas components, vector slip, Re and initial drag/net acceleration.
data/initial_interpolation_stencils.json contains the original corner values.
The exact crossing falls between PZT rows 202
and 201; both components are linearly interpolated
in height at each native column, with no adjacent-column averaging. These are
developed droplets at the interaction height, not droplets reset to their initial
near-atomiser velocity or relocated atomiser exits. Vector sample counts are not
interpreted as droplet concentration.

Reference release: droplet (u,v)=(0.0733083814,
0.964010941) m/s, sampled gas (u,v)=
(-0.325258338,0.961214737) m/s.
The plotted gas centreline is 0.964010941 m/s;
sampled-minus-plotted axial difference is -0.00279620461 m/s.
Initial net axial acceleration is -19.4884592 m/s2.
The crossing equates axial profiles, not vector velocities. In ZNMF pipeline.py,
extract_centreline finds subgrid peak positions and linearly interpolates velocities
there. Interpolating that moving-centreline curve downstream differs from bilinear
sampling at a single interpolated physical transverse position. Off-axis releases
add real transverse-profile differences. No gas value is changed to force equality.

## Actual predictions and common supported comparison

Completed cases: 7. Terminations:
{'outside_measured_jet_domain': 7}. All histories begin at the interaction plane.
The common first-ascent interval is y=10.0860744..27.3064 mm, limited by BOTH
measured PZT profile support and all computed ascending trajectories. Velocity
changes are 0.928136413..1.34178257 m/s and travel times 12.710875..13.9472779 ms.
Reference release: delta_v=1.04453009 m/s over
13.2889135 ms, ending v=2.00854103 m/s.
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
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook

report_manifest.json explicitly selects active PDF/600-dpi PNG figures and data.
source_manifest.json records original paths/configuration, transformations and hashes.
interaction_plane_audit.json records source reproduction, the computed root, and
pre-integration checks. data/ contains initial audits, full trajectories, the exact
main-comparison sources, common-height comparisons, force intervals, scalars and
local dimensionless diagnostics. validation/ contains convergence/preservation checks.
Archived H/L/C and superseded near-atomiser predictions are excluded from the notebook.
