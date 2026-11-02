# Upward co-directed prescribed-field prediction

Actual completed integrations: 7. The experiments
were separate; no combined-flow validation, collection efficiency or lung-deposition
claim is made. The obsolete horizontal side-jet study is not regenerated.

## Geometry and measured initial states

Prediction x is horizontal/rightward and y upward. The documented PZT fixed-line
reference is x=14.54665 mm, native column 42,
through the global mean-v maximum (methods/vertical_profile_selection.json in the
authoritative PZT results). It is an empirical plume reference, not a located outlet.
Nominal gas axis is aligned to it; real gas-centreline wandering is not straightened.
Gas rotation is +90 degrees: (x,y)=(x_ref-Yj, Xj+translation_y),
(u_g,v_g)=(-V_ZNMF,+U_ZNMF). No source velocity is scaled or sign-tuned.
The PZT lowest measurement plane y=0.0856 mm and first
supported gas downstream plane Xj=0.3232 mm are aligned
with axial offset 0 mm. This is an assumed coincidence
of MEASUREMENT PLANES, not physical exits. Offset stays configurable in co_directed.
All transformations apply to analysis copies. Exact cardinal rotations avoid
floating-point boundary artefacts; no measured-field extension is introduced.

Seven native PZT columns [39, 40, 41, 42, 43, 44, 45] at row 318 supply
the first common supported release plane, with both signed measured components.
All starts and source indices are in data/trajectory_summary.csv. These are
representative measured velocity samples, not tracked droplet identities or
droplet concentration. Vector sample counts are never used as number concentration.
No unmeasured travel from the atomiser surface is simulated. PZT velocities are
used once for initialisation and are never added to gas velocities or imposed again.

## Actual results

Initial sampled gas axial velocities range -0.6891127
to -0.1950244 m/s. Negative local axial gas
measurements near the nominal axis at the first plane are preserved, despite the
nominal upward jet direction. They are not masked or made positive to force coflow.
Cases ending with downward droplet velocity: 7.
Cases with any upward drag: 0; cases with any net
upward axial acceleration: 0.
Maximum heights span 0.17518969 to
0.33454286 mm. Terminations: {'outside_measured_jet_domain': 7}.

Comparisons use the SAME first-ascent height interval y=0.0856 to
0.17518969 mm, not differing endpoints. Travel times are 0.06491922 to
0.2775 ms; axial velocity changes are -0.7710215 to -0.4076291 m/s.
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
and Stk_D is zero. Stk_f=tau_Stokes*f=0.0366990043 uses documented jet actuation
f=120 Hz, not acquisition 10 Hz or atomiser 190 kHz, and
does not include 2*pi. Scalars are in data/scalar_quantities.csv, not constant-line plots.
Local gas Re_D,c=Uc*D0/nu and Sr_D,c=f*D0/Uc use D0=2 mm and positive finite
measured Uc, retaining gaps. Stk denotes Stokes number; Sr denotes Strouhal number.
No nozzle-exit Re0 or Sr0 is inferred. Supporting PZT-only Re uses an explicit
ZERO-BACKGROUND-GAS REFERENCE, not a measurement of ambient PZT gas velocity.

## Reproduce and report

From PZT_droplets workspace root:
`.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook`

Active report selection is report_manifest.json. figures/ contains seven title-free
PDF and 600-dpi PNG sets. methods_results.tex and captions.tex are report snippets.
data/ contains seven histories, trajectory summary, common-interval comparison,
force intervals, scalar table, local gas numbers, PZT-only baseline and explicitly
conditional PZT reference Re. Validation records cover common-frame rotation,
force signs, support, step convergence and unchanged source hashes. The archive
manifest is PZT_Atomiser_Analysis\Archive\20261102T103931_AEDT_superseded_side_jet\archive_manifest.json; old H/L/C identifiers and
constant-line figures are excluded. Source paths, processing configurations,
original frame definitions and applied analysis-copy transforms are in source_manifest.json.
