"""Report the measured-axis estimate without interpreting it as a remounting experiment."""


def write_report(out,cfg,metadata,stations,figures,checks):
    fit=metadata['primary_fit'];sens=metadata['fit_sensitivity'];origin=metadata['reference_point_lab_mm']
    angle=fit['angle_right_of_lab_vertical_deg'];angles=[r['angle_right_of_lab_vertical_deg'] for r in sens]
    metrics=[r for r in stations if r['width_asymmetry'] is not None]
    text=f'''# Rigid plume-aligned measured PZT field

Interaction-prediction work is paused. This package does not run a prediction,
select an interaction station or change authoritative experimental data.

## Orientation evidence and fitted axis

RAW_data/Batch4_Cam1.h5 is the available PIV input; no linked original PIV image
sequence or acquisition-angle record was found. Nearby Atomiser images/VMARF
TIFFs identify an EVOS M5000 microscope, 2048x1536 pixels, about 0.619 um/pixel,
and single microscope frames, rather than this 5.35 um/pixel PIV field. The
assembly illustration is not an acquisition-angle measurement. See
validation/image_inventory.json and its inspected preview. Existing processing
uses us/vs as adopted horizontal/vertical displacements, retains their signs and
applies 0.214 m/s per pixel/pair once before statistics. This analysis reads the
saved SI mean/covariance fields directly and applies no calibration again.

One hydrodynamic axis is estimated from the mean laboratory vertical component.
At each native laboratory y plane a three-node concave quadratic locates the
subgrid maximum. The connected half-maximum lobe must contain at least
{cfg['minimum_connected_halfmax_native_nodes']} native nodes, its peak must exceed
{cfg['minimum_positive_peak_m_s']} m/s, and the required paired values must be finite.
No smoothing of the field, profile or extracted centres is applied. Fit weights
are equal per native y plane; sample counts are not treated as independent droplet counts.

The primary laboratory interval requested is {cfg['fit_interval_lab_y_mm']} mm;
actual accepted y is {fit['actual_lower_y_mm']:.8g} to {fit['actual_upper_y_mm']:.8g} mm,
with {fit['centre_count']} centres. Half-maximum lobes span
{fit['halfmax_nodes_min']}--{fit['halfmax_nodes_max']} native nodes. OLS fits
x_c = a + b*y, a={fit['intercept_x_at_lab_y0_mm']:.10g} mm,
b={fit['slope_dx_dy']:.10g}. The axis leans **{angle:.7g} degrees right of laboratory
vertical**. Residual RMS is {fit['residual_RMS_x_mm']:.7g} mm and maximum absolute
residual {fit['residual_max_abs_x_mm']:.7g} mm over that interval. These are fit
residuals of subgrid-localised flow peaks, not instrumental accuracy or independent
sample confidence limits. Native spacings remain 0.3424 mm horizontally and
0.0856 mm vertically; interrogation resolution is not inferred from them.

Other requested intervals give angles {min(angles):.7g} to {max(angles):.7g} degrees.
The full table is data/fit_interval_sensitivity.csv. A connected-half-maximum-lobe
velocity-weighted centroid provides a separate estimator sensitivity check in
transformation.json; it does not replace or recenter the data. Neither curvature
nor asymmetry is assumed to come entirely from mounting tilt. The angle is a
measured-flow straight-axis estimate, not an independently verified mechanical angle.

The secondary estimator gives {metadata['secondary_centre_estimator_fit']['angle_right_of_lab_vertical_deg']:.7g}
degrees over the primary interval, providing a separate method sensitivity check.

## One rigid frame transformation

Reference point in laboratory coordinates: (x,y)=({origin[0]:.10g}, {origin[1]:.10g}) mm.
It is the fitted axis at the lowest accepted native plane in the fitting interval,
inside measured support. It is NOT the atomiser surface or physical outlet.
For theta={angle:.9g} degrees, Q=[[cos(theta),-sin(theta)],
[sin(theta),cos(theta)]]. Coordinates [n_p,s_p]=Q*([x,y]-reference);
velocities [u_p,v_p]=Q*[u,v]. n_p is transverse/right of the fitted axis and
s_p increases along the plume. The s_p axis is PLUME-ALIGNED, not laboratory
vertical; points below the reference legitimately have negative s_p.
The same determinant-plus-one orthonormal matrix acts on positions and velocity
components, preserving magnitudes. The inverse is Q transpose plus translation.
Laboratory gravity [0,-9.80665] m/s2 becomes
{metadata['gravity_plume_m_s2']} m/s2 in [n_p,s_p] components, saved for future models.
No interaction model is run here. Re-expressing one measured field in this frame
is not evidence that remounting the physical atomiser gives identical flow.

Every rotated native node is retained without resampling for contour maps.
The field is not bent; rows are not shifted independently; neither side is mirrored,
averaged with the other or forced symmetric. Transverse velocity is not forced zero.
The original centre positions retain nonzero n_p residuals and their curvature.

## Statistics, axis profile and transverse sections

The complete saved velocity covariance tensor is transformed as C_p=Q*C_lab*Q^T.
Aligned RMS values are square roots of its diagonal, not RMS values rotated as
vectors. Covariance/correlation and unchanged magnitude scalars are saved in
data/statistics_plume_aligned.npz, together with native indices/support and N.
The original finite-pair population and N remain unchanged.

The main axial profile follows the SINGLE fitted straight axis n_p=0 over measured
support, sampled at the native laboratory y planes using finite-corner bilinear
interpolation in the original rectangular grid. This samples a different spatial
line from the original fixed laboratory column. The original CSV and NPZ are copied
byte-for-byte under reference/. Their original axes/units remain intact. A separate
comparison re-expresses original fixed-column points by s_p and shows both original
laboratory and projected axial components, separating component rotation from line
selection. No missing values are filled and no profile is extrapolated.

Transverse sections are at fixed plume-aligned stations s_p={cfg['transverse_stations_plume_s_mm']} mm.
Their original laboratory positions and support flags are exported. These cuts
are perpendicular to the one fitted axis, not individually recentered profiles.
The half-width metrics use each cut's actual peak and its first supported
half-maximum crossings; neither side is altered to form them. Width asymmetry is
(right_half_width-left_half_width)/(right_half_width+left_half_width).
Transverse-profile plots show n_p between -7 and +7 mm for plume detail; full
supported cuts are retained in CSV, and all native data remain in the contour arrays.

'''
    text+='| s_p (mm) | peak n_p (mm) | left half-width (mm) | right half-width (mm) | width asymmetry |\n|---|---|---|---|---|\n'
    for r in stations:
        def fmt(v):return 'missing' if v is None else f'{v:.6g}'
        text+='| '+' | '.join(fmt(r[k]) for k in ['s_p_mm','peak_n_mm','left_half_width_mm','right_half_width_mm','width_asymmetry'])+' |\n'
    text+='''
These differences and the retained nonzero transverse component are descriptive
evidence of residual asymmetry in the measured representation, not a quantified
instrumental-significance test. No mounting-only explanation, symmetry correction
or uncertainty smaller than a grid cell is claimed from subgrid localisation.

## Reproduce and inspect

From the PZT_droplets workspace root:

```powershell
.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_plume_aligned_notebook
```

Reusable code: analysis/plume_alignment.py and plume_alignment_report.py.
Configuration: config/plume_alignment.json. The active notebook displays this
alignment section while archived/preserved interaction displays are paused.
report_manifest.json selects only these measurement figures; transformation.json
contains Q, inverse, origin, gravity, fitting interval and sensitivity. Original
source paths and calibration provenance are in source_manifest.json. Numerical
checks verify inverse positions/vectors, velocity-magnitude preservation,
covariance inverse/eigenvalues, finite support, N and protected file hashes.
'''
    (out/'README.md').write_text(text,encoding='utf-8')
    tex=rf'''\section{{Rigid plume-aligned representation of the measured PZT field}}
No linked acquisition image establishes the mounting orientation. Nearby EVOS
microscope images and an assembly illustration are not calibrated orientation
evidence for this PIV sequence. A straight hydrodynamic axis is therefore fitted
to subgrid maxima of the laboratory mean vertical component over
$y={fit['actual_lower_y_mm']:.7g}$--${fit['actual_upper_y_mm']:.7g}\,$mm
({fit['centre_count']} native planes; {fit['halfmax_nodes_min']}--{fit['halfmax_nodes_max']}
nodes across each connected half-maximum lobe). No field smoothing is applied.
The fitted axis leans ${angle:.7g}^\circ$ right of laboratory vertical;
RMS and maximum absolute horizontal residuals are {fit['residual_RMS_x_mm']:.6g}
and {fit['residual_max_abs_x_mm']:.6g} mm. Fits over alternative intervals span
{min(angles):.6g}--{max(angles):.6g} degrees. This sensitivity is not a confidence
interval, and the estimate does not establish that curvature or asymmetry is caused
solely by mounting tilt.

The reference point $(x_r,y_r)=({origin[0]:.8g},{origin[1]:.8g})\,$mm lies on
the fitted axis within the fitting interval; it is not the atomiser surface.
With $\theta={angle:.8g}^\circ$,
\begin{{align}}
Q&=\begin{{pmatrix}}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{{pmatrix}},\\
\begin{{pmatrix}}n_p\\s_p\end{{pmatrix}}&=Q\left(\begin{{pmatrix}}x\\y\end{{pmatrix}}-
\begin{{pmatrix}}x_r\\y_r\end{{pmatrix}}\right),\qquad
\begin{{pmatrix}}u_p\\v_p\end{{pmatrix}}=Q\begin{{pmatrix}}u\\v\end{{pmatrix}}.
\end{{align}}
$s_p$ is plume-aligned axial distance, not laboratory vertical. The rotated native
grid is plotted directly without resampling. One rigid transform retains curvature,
asymmetry and transverse motion; no row-dependent recentering or symmetry operation
is applied. Velocity magnitudes are preserved. Covariance is transformed as
$C_p=QCQ^T$, and RMS follows its diagonal. Gravity is saved as $Q(0,-g)^T$ for
future modelling. This coordinate change does not simulate physical remounting.

The fitted-axis profile samples $n_p=0$, a different spatial line from the original
fixed-column profile, which is preserved separately. Axis and transverse sections
use finite-corner bilinear interpolation only within original measured support;
missing values and boundaries are not extrapolated. Transverse stations have
fixed $s_p$, rather than independent local centres. Half-width asymmetry and
nonzero transverse velocities remain descriptive measured features; subgrid fit
residuals are not interpreted as instrumental accuracy. Original SI calibration,
sample populations and saved experimental results are unchanged. Interaction
predictions are paused and no interaction station is selected here.
'''
    (out/'methods.tex').write_text(tex,encoding='utf-8')
    captions={
        'original_field_axis_fit':'Original laboratory mean vertical velocity, extracted subgrid peak centres and the single fitted straight axis. Black centres lie within the fitting interval and grey centres outside. The open marker is the documented analysis reference point, not the atomiser surface. The field is unsmoothed.',
        'plume_aligned_mean_axial_velocity':'Mean plume-aligned axial velocity on the rotated native grid. n_p is transverse to the one fitted axis and s_p increases along the plume, not laboratory vertical. Signed values, native support and asymmetry are retained without resampling or row-dependent shifts.',
        'plume_aligned_mean_transverse_velocity':'Mean plume-aligned transverse component after the same orthonormal rotation of both velocity components. No zero-transverse constraint or side averaging is imposed; remaining transverse motion and asymmetry are measured features of this representation.',
        'fitted_axis_axial_velocity_profile':'Mean plume-aligned axial velocity sampled along the fixed fitted straight axis n_p=0. This samples a different spatial line than the original laboratory fixed-column profile. Finite-corner interpolation uses only measured support, preserving gaps without extrapolation.',
        'transverse_axial_velocity_profiles':'Plume-aligned axial velocity at fixed documented s_p stations on cuts perpendicular to the single fitted axis. The displayed transverse window is minus/plus 7 mm; full supported cuts are exported. Neither side is mirrored, averaged or recentered.',
        'transverse_transverse_velocity_profiles':'Plume-aligned transverse velocity on the same fixed-s_p sections, retaining signed values and asymmetry. No transverse velocity is forced to zero.',
        'original_vs_plume_aligned':'Original laboratory vertical component and plume-aligned axial component, each on its native physical coordinates. Both panels use the same linear symmetric colour limits. A rigid coordinate/component transformation changes representation, not the original measurements or the physical mounting.',
        'different_profile_lines_comparison':'Original fixed-column laboratory component, its plume-aligned component at the same physical points, and the aligned component along the fitted axis. Fixed-column points are re-expressed by s_p only for comparison. The curves separate component rotation from the effect of sampling a different spatial line.',
        'straight_axis_residuals':'Extracted centre offsets from the one fitted straight axis versus laboratory y. The shaded interval identifies the primary fit region. Structured residuals and interval sensitivity can reflect genuine curvature, asymmetry and estimator behaviour, not solely mounting tilt.'}
    (out/'captions.tex').write_text(''.join('% '+stem+'\n'+r'\caption{'+captions[stem]+'}\n\n' for stem in figures),encoding='utf-8')
