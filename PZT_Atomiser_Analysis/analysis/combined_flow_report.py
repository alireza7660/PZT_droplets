"""Methods and report selection for separate experiments and a proposed side jet."""


def write_report(out,cfg,sources,status,figures):
    missing='\n'.join('- '+s for s in status['missing_inputs']) or '- None; explicit configured cases have been evaluated.'
    diameter_note=cfg.get('diameter_source_or_assumption') or 'No diameter has been configured.'
    text=f'''# Proposed combined side-jet / droplet-plume workflow

## What was measured

PZT-only: upward droplet plume, atomiser at 190 kHz and 9 Vpp, ZNMF jet absent.
The PIV values are droplet-phase velocities. Surrounding-air velocity was not
measured; this experiment is not labelled stationary air. ZNMF-only: gas velocities
measured without PZT droplets. The proposed combined experiment has not been conducted.

The active status is **{status['status']}**. Baseline results are complete.
No arbitrary diameter or placement cases are populated. No combined-flow prediction,
entrainment, capture or experimental model validation is claimed while inputs are missing.

## Specific missing prediction inputs

{missing}

These are the only presently missing experimental/design inputs for the implemented
first-order, two-dimensional model. A validated size distribution is welcome, but
an explicitly declared diameter study is also supported. No diameter is inferred
from pixel calibration, interrogation-window size, image-spot width or nozzle diameter.
No numerical droplet prediction is substituted for these missing inputs.

## Configured nominal droplet diameter

{diameter_note}

The adopted 10 micrometre value, not the unrounded 9.45 micrometre estimate, is
the nominal input to drag, particle Re, Stokes-limit response time and local Stokes
number. The aperture-to-jet equality and classical cylindrical-jet breakup are
modelling assumptions, not independently verified atomiser sizing. This does not
establish a measured droplet size, monodispersity or a validated size distribution.
`data/droplet_model_parameters.csv` records the configured diameter and response
time. Size and provenance remain configurable for measured sizing or a declared
later parameter study. Prediction Re and St are not reported until gas slip and
placement/release inputs can be evaluated.

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
The two native-column profiles use the existing fixed x={status['fixed_profile_maximum']['x_mm']:.8g} mm
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
.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_combined_flow_notebook
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
old figures. The timestamped archive manifest is `{cfg['archive_manifest']}`.
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
'''
    if cfg.get('design_study',{}).get('enabled') and status.get('prediction_case_count',0)>0:
        text=text.replace('No arbitrary diameter or placement cases are populated. No combined-flow prediction,\nentrainment, capture or experimental model validation is claimed while inputs are missing.',
            f"The prospective design study integrated {status['prediction_case_count']} conditional particle cases in {status.get('proposed_design_placement_count',0)} proposed arrangements. Placement is a DESIGN VARIABLE, not measured geometry. See Design_Study/placement_comparison.csv and case_comparison.csv. No capture efficiency, optimal arrangement or experimental validation is claimed.")
        text=text.replace('No numerical droplet prediction is substituted for these missing inputs.',
            'Fixed measured/prototype geometry is not required for the proposed placement study; it is selected and documented from measured coordinate support.')
        text=text.replace('## Specific missing prediction inputs','## Fixed-geometry inputs (not prerequisites for the prospective design study)')
        text+='\n## Active prospective outputs\n\nDesign_Study/manifest.json lists proposed trajectories, velocity/displacement/Re/Stokes figures and comparison tables. Dimensionless_Numbers/manifest.json lists independent PZT-reference and ZNMF local-number CSV/PDF/PNG outputs. Dimensionless calculations can also run without any geometry using `python -m PZT_Atomiser_Analysis.analysis.prospective_studies --dimensionless-only`.\n'
    (out/'README.md').write_text(text,encoding='utf-8')
    tex=r'''\section{Prediction workflow for the proposed side-jet arrangement}
The PZT-only experiment measured droplet-phase velocities in an upward plume
without the ZNMF jet. Its surrounding-air velocity was not measured. The separate
ZNMF-only experiment measured a gas jet without PZT droplets. The proposed combined
experiment has not been conducted. Consequently, the experimental fields are
kept separate and predictions are not described as validated combined measurements.

Final PZT velocities retain the $5.35\,\mu\mathrm{m/pixel}$ calibration,
$25\,\mu\mathrm{s}$ pulse separation and $0.214\,\mathrm{m\,s^{-1}}$ per pixel
displacement conversion already applied in processing. The atomiser operates at
190 kHz and 9 Vpp. Signed mean fields and existing sample selection are unchanged.
The fixed vertical-line profiles are ensemble means, not particle trajectories
or tracked plume centrelines. No ambient gas velocity or measured droplet size
is inferred from them.

The common frame retains the PZT image-reference origin, $x$ rightward and $y$
upward. The ZNMF registered coordinates are dimensionalised once using $D_0=2\,$mm
and translated to its documented geometric nozzle centre. Configured placement
then applies $\boldsymbol q=R\boldsymbol q_j+\boldsymbol q_N$ and
$\boldsymbol w_g=R(U,V)$ to positions and components consistently; independent
experimental origins are not equated.

When diameter and placement/release inputs are explicitly provided, PZT velocities
initialise representative droplets once. They are not subsequently added to the
gas or repeatedly imposed on the particles. The separately measured two-component
ZNMF mean alone is the prescribed gas field:
\begin{align}
\dot{\boldsymbol q}_d&=\boldsymbol w_d,\\
\dot{\boldsymbol w}_d&=-g(1-\rho_a/\rho_w)\boldsymbol e_y+
\frac{f(Re_d)}{\tau_{p,\mathrm{Stokes}}}(\boldsymbol w_g-\boldsymbol w_d),\\
Re_d&=\frac{\rho_a d\lVert\boldsymbol w_g-\boldsymbol w_d\rVert}{\mu_a},
\qquad \tau_{p,\mathrm{Stokes}}=\frac{\rho_wd^2}{18\mu_a}.
\end{align}
Here $f=1$ for $Re_d\leq0.01$, otherwise $f=1+0.15Re_d^{0.687}$, restricted
to $Re_d<800$. The force is finite at zero slip. Gas evaluation is confined to
four-finite-vector measured cells, without extrapolation; predictions terminate
at unsupported boundaries/gaps, a drag-law limit or the configured time horizon.
Gravity acts downward in the physical frame; jet rotation does not rotate gravity.

The model prescribes a gas-only mean and neglects plume-induced gas changes and
jet unsteadiness. Constant spherical diameter, no evaporation and one-way coupling
are assumptions. At 20 degrees C, dry-air properties approximate unknown humidity;
101.325 kPa and equal water/air temperature are assumed. Property sources and
values are recorded separately. No dominant entrainment or transport mechanism
is inferred without the required inputs or combined measurements.

The Stokes-limit response time is distinct from the finite-Re drag/slip coefficient
timescale $\tau_{\mathrm{drag,coef}}=\tau_{p,\mathrm{Stokes}}/f(Re_d)$.
An explicitly defined local comparison may use $\tau_{flow}=D_0/\lVert\boldsymbol w_g\rVert$
and $St_{D,local}=\tau_{p,\mathrm{Stokes}}\lVert\boldsymbol w_g\rVert/D_0$.
This is not a nozzle-exit or atomiser-frequency response comparison and cannot
alone predict transport, capture or successful entrainment.
'''
    tex+=r'''
The nominal mesh aperture diameter is $5\,\mu\mathrm{m}$. Assuming that the
emerging cylindrical jet has the same diameter, the classical Rayleigh estimate
$d_{drop}\simeq1.89d_{jet}$ gives $9.45\,\mu\mathrm{m}$. A nominal
$d=10\,\mu\mathrm{m}$ is adopted for the initial case. This is a model-based
estimate, not a measured droplet size or validated size distribution. Neither
the aperture-to-jet equality nor Rayleigh breakup applicability is independently
validated for this driven mesh atomiser. The adopted value, rather than the
unrounded estimate, is used consistently for drag, Reynolds number,
$\tau_{p,\mathrm{Stokes}}$ and Stokes number; diameter remains configurable
for replacement by measured sizing. The configured response time is exported
even when geometry is not yet sufficient for transport prediction.
'''
    if status['missing_inputs']:
        tex+='\nMeasured baselines, the nominal response-time table and the reusable input-gated workflow are reported. Numerical combined-flow predictions await the remaining inputs listed in prediction_status.json, currently side-jet nozzle placement/direction and PZT release locations. No placement is invented.\n'
    (out/'methods.tex').write_text(tex,encoding='utf-8')
    captions=''
    for item in figures:captions+='% '+item['stem']+'\n'+r'\caption{'+item['caption'].replace('_',r'\_')+'}\n\n'
    (out/'captions.tex').write_text(captions,encoding='utf-8')
