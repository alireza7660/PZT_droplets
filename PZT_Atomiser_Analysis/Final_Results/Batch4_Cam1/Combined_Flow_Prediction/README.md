# Proposed combined side-jet / droplet-plume workflow

## What was measured

PZT-only: upward droplet plume, atomiser at 190 kHz and 9 Vpp, ZNMF jet absent.
The PIV values are droplet-phase velocities. Surrounding-air velocity was not
measured; this experiment is not labelled stationary air. ZNMF-only: gas velocities
measured without PZT droplets. The proposed combined experiment has not been conducted.

The active status is **prospective_design_study_complete**. Baseline results are complete.
The prospective design study integrated 27 conditional particle cases in 9 proposed arrangements. Placement is a DESIGN VARIABLE, not measured geometry. See Design_Study/placement_comparison.csv and case_comparison.csv. No capture efficiency, optimal arrangement or experimental validation is claimed.

## Fixed-geometry inputs (not prerequisites for the prospective design study)

- None; explicit configured cases have been evaluated.

The completed prospective study does not require measured combined geometry. Proposed placements and supported release nodes are specified numerically in its comparison tables. A validated size distribution is welcome, but
an explicitly declared diameter study is also supported. No diameter is inferred
from pixel calibration, interrogation-window size, image-spot width or nozzle diameter.
Fixed measured/prototype geometry is not required for the proposed placement study; it is selected and documented from measured coordinate support.

## Configured nominal droplet diameter

User-selected nominal model estimate: 5 um nominal mesh aperture, emerging jet diameter assumed equal to aperture; Rayleigh d_drop approximately 1.89*d_jet = 9.45 um, rounded to 10 um. Not measured sizing or a validated size distribution.

The adopted 10 micrometre value, not the unrounded 9.45 micrometre estimate, is
the nominal input to drag, particle Re, Stokes-limit response time and local Stokes
number. The aperture-to-jet equality and classical cylindrical-jet breakup are
modelling assumptions, not independently verified atomiser sizing. This does not
establish a measured droplet size, monodispersity or a validated size distribution.
`data/droplet_model_parameters.csv` records the configured diameter and response
time. Size and provenance remain configurable for measured sizing or a declared
later parameter study. Trajectory slip-based Re and local Stokes diagnostics are reported in Design_Study. Geometry-independent reference and local-centreline numbers are reported in Dimensionless_Numbers.

## Authoritative data and physical frames

`source_manifest.json` gives full source paths, hashes, coordinate formulas and
mask policies. PZT mean_u/mean_v and saved frame u/v already have units m/s;
they are read unchanged. PZT calibration remains 5.35 micrometres/pixel,
25 microseconds and 0.214 m/s per pixel displacement, applied once in the original
processing. Existing signed values, sample counts and sample selection are preserved.
PZT x is horizontal/rightward and y vertically upward in the exported image-reference
frame. That origin is not an established atomiser surface.

The full ZNMF-only reconstructed field is used, not the old cropped centreline CSV.
ZNMF U,V are already m/s, including the upstream processing sign convention.
No additional calibration or sign reversal is applied. Its stored coordinates are
dimensionless: Xj=x*D0 and Yj=(y-6.868)*D0, with D0=0.002 m and the documented
geometric nozzle centre at registered y/D0=6.868. These are local nozzle-centred
coordinates; source y_plot uses a different velocity-peak reference and is not used.
Only coordinate nondimensionalisation is reversed, not velocity calibration.
The source gas reconstruction retains its existing FOV registration/blending.

The common frame retains the PZT image-reference axes and origin. An explicitly
configured nozzle position and rotation transform both gas positions and components:
q_common=R(theta)*q_jet+q_nozzle; w_common=R(theta)*(U,V).
Translation does not act on velocity. Nozzle and plume origins are never silently
identified. A horizontal side jet may enter from either side (rotation 0 or 180 degrees
in the PZT frame); the actual side, nozzle position and release location are user inputs.
No gas/PZT overlay is made before those inputs are supplied.

## Measured baselines

`pzt_only_mean_u` and `pzt_only_mean_v` show the original signed mean fields.
The two native-column profiles use the existing fixed x=14.54665 mm
through the global mean-v maximum. All rows, signs and gaps are retained without
smoothing, adjacent-column averaging or extrapolation. They are spatial ensemble
means, not tracked individual droplets or a tracked plume centreline. A lateral
plume shift can change a fixed-line profile with height. The separate ZNMF-only
component maps show input support in their own local frame and imply no placement.

## Prescribed-field prediction when inputs are supplied

For each explicitly supplied diameter and native PZT release node, measured
droplet components initialise position and velocity. `velocity_source` is
`ensemble_mean` or `frame:<saved key>`; these initialise representative model
particles and do not establish individual tracking. The PZT-only field is never
imposed again after release. Neither particle velocity nor the PZT field is added
to the gas velocity. The two-component measured ZNMF mean field alone supplies
gas velocity at each evolving particle position.

Particle acceleration is -g*(1-rho_air/rho_water)*e_y +
(w_gas-w_droplet)*f(Re)/tau_Stokes, with Re=rho_air*d*norm(vector slip)/mu_air,
tau_Stokes=rho_water*d^2/(18*mu_air), f=1 for Re<=0.01 and
f=1+0.15*Re^0.687 otherwise. The Schiller-Naumann implementation is restricted
to Re<800, using the cited COMSOL range. Zero slip is finite. All integration is
in SI units; gravity remains downward in the common frame, independent of jet rotation.

Gas evaluation uses bilinear interpolation only within source-grid cells with
four finite paired components. It does not smooth the experimental baselines or
extrapolate gas velocity. RK4 stages and traversed segments are checked; predictions
stop at the last supported state at a source boundary or invalid cell, a drag-law
limit or the configured time horizon. An upward turning point is not a stop condition
for side-jet transport: later deflection/settling can occur while measured support remains.
Half-step comparisons are recorded separately in validation/prediction_convergence.json.

This is a first-order one-way prescribed-mean-field model. It neglects changes
to the gas caused by the plume, jet unsteadiness, two-way coupling, evaporation,
diameter changes, collisions, out-of-plane motion, internal circulation and deformation.
Air is at the recorded 20 C; 101325 Pa and equal water/air temperature are assumptions.
Humidity is unavailable; dry-air properties are an approximation. Property values
and authoritative source URLs are recorded in fluid_properties.json.

## Reynolds, response time and Stokes definitions

Predicted Re uses the evolving in-plane vector gas-droplet slip and each explicit
diameter. No measured PZT-only Re is claimed because ambient gas is unknown and the configured size is a model estimate.
tau_Stokes is the low-Re response-time baseline. The exported finite-Re quantity
tau_drag_coefficient=tau_Stokes/f(Re) is an instantaneous drag/slip coefficient
timescale, not an exact finite-Re relaxation time. The optional local response
comparison uses tau_flow=D0/norm(local prescribed mean gas velocity) and
St_D_local=tau_Stokes*norm(gas velocity)/D0. At zero gas magnitude tau_flow is
infinite and St_D_local is zero. This chosen local nozzle-length timescale is not
a nozzle-exit or 190 kHz forcing timescale. Stokes number alone does not predict
droplet transport, capture or successful entrainment.

## Reproduce and configure

From the PZT_droplets workspace root:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook
```

This runs the reusable workflow, refreshes only the notebook's prediction section,
and validates outputs. It does not rerun the experimental accumulation. Configuration:
`PZT_Atomiser_Analysis/config/combined_flow_prediction.json`.
Set `droplet_diameters_um` and `diameter_source_or_assumption`; set geometry's
nozzle position [x_mm,y_mm], rotation and placement_source; set release_locations
to objects containing row, column and velocity_source. No example numerical cases
are pre-populated. Invalid or unsupported release locations are reported, not moved
to make a result. The common frame and field masks remain explicit.

The underlying command-line engine is analysis/combined_flow_prediction.py;
measured-gas interpolation and integration are in analysis/prescribed_gas.py;
force and property functions are in analysis/droplet_dynamics.py.

## Report and archive

Use only the explicit report_manifest.json entries; do not recursively discover
old figures. The timestamped archive manifest is `PZT_Atomiser_Analysis/Archive/20261102T095744_AEDT_superseded_interaction/archive_manifest.json`.
It records every moved theory/jet-Stokes/axial-overlay file with original relative
path and SHA-256, source/configuration snapshots and exact replacement mapping.
Original inputs, authoritative PZT/ZNMF results, unrelated analyses and Git history
are excluded from the move. The old arbitrary-size and same-axis comparisons have
no new numerical counterpart until required prediction inputs are specified.

Superseding outputs are the measured baseline figures/CSV/NPZ, authoritative
source manifest, input-gated configuration and, only when configured, integrated
trajectory/component/slip/Re/St CSV and figures. `prediction_status.json` is the
authoritative statement of what was generated. Detailed synthetic solver tests
are validation checks, not experimental parameter-study results.

## Active prospective outputs

Design_Study/manifest.json lists proposed trajectories, velocity/displacement/Re/Stokes figures and comparison tables. Dimensionless_Numbers/manifest.json lists independent PZT-reference and ZNMF local-number CSV/PDF/PNG outputs. Dimensionless calculations can also run without any geometry using `python -m PZT_Atomiser_Analysis.analysis.prospective_studies --dimensionless-only`.
