# Rigid plume-aligned measured PZT field

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
7 native nodes, its peak must exceed
0.5 m/s, and the required paired values must be finite.
No smoothing of the field, profile or extracted centres is applied. Fit weights
are equal per native y plane; sample counts are not treated as independent droplet counts.

The primary laboratory interval requested is [5.0, 22.0] mm;
actual accepted y is 5.0504 to 21.9992 mm,
with 199 centres. Half-maximum lobes span
9--16 native nodes. OLS fits
x_c = a + b*y, a=14.50303181 mm,
b=0.09204581548. The axis leans **5.259018 degrees right of laboratory
vertical**. Residual RMS is 0.03261133 mm and maximum absolute
residual 0.07902312 mm over that interval. These are fit
residuals of subgrid-localised flow peaks, not instrumental accuracy or independent
sample confidence limits. Native spacings remain 0.3424 mm horizontally and
0.0856 mm vertically; interrogation resolution is not inferred from them.

Other requested intervals give angles 4.864496 to 5.431601 degrees.
The full table is data/fit_interval_sensitivity.csv. A connected-half-maximum-lobe
velocity-weighted centroid provides a separate estimator sensitivity check in
transformation.json; it does not replace or recenter the data. Neither curvature
nor asymmetry is assumed to come entirely from mounting tilt. The angle is a
measured-flow straight-axis estimate, not an independently verified mechanical angle.

The secondary estimator gives 5.232024
degrees over the primary interval, providing a separate method sensitivity check.

## One rigid frame transformation

Reference point in laboratory coordinates: (x,y)=(14.9679, 5.0504) mm.
It is the fitted axis at the lowest accepted native plane in the fitting interval,
inside measured support. It is NOT the atomiser surface or physical outlet.
For theta=5.25901793 degrees, Q=[[cos(theta),-sin(theta)],
[sin(theta),cos(theta)]]. Coordinates [n_p,s_p]=Q*([x,y]-reference);
velocities [u_p,v_p]=Q*[u,v]. n_p is transverse/right of the fitted axis and
s_p increases along the plume. The s_p axis is PLUME-ALIGNED, not laboratory
vertical; points below the reference legitimately have negative s_p.
The same determinant-plus-one orthonormal matrix acts on positions and velocity
components, preserving magnitudes. The inverse is Q transpose plus translation.
Laboratory gravity [0,-9.80665] m/s2 becomes
[0.898861356725992, -9.765369039820497] m/s2 in [n_p,s_p] components, saved for future models.
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

Transverse sections are at fixed plume-aligned stations s_p=[0.0, 5.0, 10.0, 15.0, 20.0] mm.
Their original laboratory positions and support flags are exported. These cuts
are perpendicular to the one fitted axis, not individually recentered profiles.
The half-width metrics use each cut's actual peak and its first supported
half-maximum crossings; neither side is altered to form them. Width asymmetry is
(right_half_width-left_half_width)/(right_half_width+left_half_width).
Transverse-profile plots show n_p between -7 and +7 mm for plume detail; full
supported cuts are retained in CSV, and all native data remain in the contour arrays.

| s_p (mm) | peak n_p (mm) | left half-width (mm) | right half-width (mm) | width asymmetry |
|---|---|---|---|---|
| 0 | 0 | 1.84776 | 1.29018 | -0.177691 |
| 5 | 0 | 2.26695 | 1.60477 | -0.171029 |
| 10 | 0 | 2.66695 | 1.92939 | -0.160467 |
| 15 | 0 | 3.02082 | 2.28701 | -0.138251 |
| 20 | 0 | 3.37218 | 2.72912 | -0.105398 |

These differences and the retained nonzero transverse component are descriptive
evidence of residual asymmetry in the measured representation, not a quantified
instrumental-significance test. No mounting-only explanation, symmetry correction
or uncertainty smaller than a grid cell is claimed from subgrid localisation.

## Reproduce and inspect

From the PZT_droplets workspace root:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_plume_aligned_notebook
```

Reusable code: analysis/plume_alignment.py and plume_alignment_report.py.
Configuration: config/plume_alignment.json. The active notebook displays this
alignment section while archived/preserved interaction displays are paused.
report_manifest.json selects only these measurement figures; transformation.json
contains Q, inverse, origin, gravity, fitting interval and sensitivity. Original
source paths and calibration provenance are in source_manifest.json. Numerical
checks verify inverse positions/vectors, velocity-magnitude preservation,
covariance inverse/eigenvalues, finite support, N and protected file hashes.
