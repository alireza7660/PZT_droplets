# PZT spatial maps

The current output is **uncalibrated diagnostic analysis**, not a physically
calibrated spray measurement. Only `RAW_data/Batch4_Cam1.h5` was supplied.
No original data are edited.

From the workspace root, run:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.pipeline --config PZT_Atomiser_Analysis/config/analysis.json
```

Open `notebooks/Spatial_maps.ipynb` using `.venv/Scripts/python.exe` as the
kernel. Its calculation and plotting cells call reusable module functions.
`analysis/plotting.py` exposes `contour_map`, `quiver_map`, and `shared_limits`.
Every exported figure has PDF, editable SVG, 600-dpi PNG and numerical NPZ.
Full spatial statistics are also in `results/Batch4_Cam1_spatial_statistics.csv`.

## Input and limitations

There are 1757 numerically named frame groups, each containing 319 x 80 arrays:
`u`, `v`, `us`, `vs`, `flag`, and `omega`. Numeric frame keys are sorted
numerically; acquisition timing and temporal order are not verified.
Native x runs from 31 to 5087 in steps of 64; native y runs from 5104 to 16
in steps of -16. This is vector spacing in unspecified native coordinate
units, not a verified interrogation-window size or spatial resolution.
No array has unit attributes. The scalar values scale=1, dt=1, calib=0,
Lc=1 and Uc=1 do not establish dimensional calibration. The full inventory
is in `results/input_inventory.json`.

Flag codes 0 and 1 occur but their meanings and replacement policy are absent.
The analysis does not infer these meanings. All paired u/v values are finite;
N=1757 at every grid location under the provisional finite-pair rule. This is
not a quality-accepted sample count. Flag-frequency plots report codes only.
There is no identifiable correlation-quality metric. `us/vs` are distinct
arrays of undocumented meaning. Their all-frame means match stored U/V
exactly, whereas the u/v means differ by up to 0.158 and 0.517 native units.
The present figures use u/v and recompute every statistic, avoiding an
undocumented substitution of us/vs. This does not prove u/v were never
processed or replaced upstream.

The user confirmed spatial calibration of 5.05 micrometres/pixel and pulse
separation of 25 microseconds. The velocity factor would be 0.202 m/s per pixel
displacement if the arrays are pixels per pair. This conversion is recorded but
not applied, since the exporting code is absent and the arrays may already be
calibrated. Likewise, coordinate conversion to mm awaits confirmation that x/y
are pixel coordinates. Confirmed calibration is not a verified storage convention.

The user also confirmed upward spray direction, upward-positive axial U,
rightward-positive transverse V, axial-position symbol x and transverse-position
symbol y. These describe the physical convention, not the unverified stored
u/v basis. Native y decreases with row index; no additional reflection is
assumed. No physical outlet position is supplied. Once mapping is verified,
figures place transverse y horizontally and axial x vertically.

Missing information: stored coordinate and velocity units/convention,
acquisition timing, verified component orientation/signs and
physical origin, flag definitions and processing software/source, relationship
between u/v and us/vs, interpolation/replacement policy, tracked phase and
operating conditions. Physical axial/transverse labels, physical-unit figures,
quality-rejected masks and report-ready physical interpretation require this
information. No forcing, calibration, nozzle or phase assumptions are inherited
from another experiment. `Report_Files/Images` remains empty until calibration
and validity are established; the generated diagnostics are clearly labelled.

## Definitions and rendering

Paired finite u/v values are used for both components. Means, RMS fluctuations
and signed covariance use population denominator N. Correlation is undefined
at zero variance. Statistics require at least two samples. The two speed
statistics are calculated separately: hypot(mean(u),mean(v)) and
mean(hypot(u,v)). No stored moments are trusted without recomputation.

Native coordinate values and equal spatial aspect are retained. No origin
shift, smoothing, extrapolation or gap filling is applied. Filled contours
interpolate contour boundaries within valid native quadrilaterals, as intrinsic
to contourf; corner_mask=False prevents contours crossing quads touching missing
nodes. This does not increase measurement resolution. No cosmetic filtering is
applied. Nonfinite paired vectors and undefined statistics are masked.

Full-range common limits are calculated before drawing the same quantity for
all configured conditions. Signed quantities use symmetric diverging scales;
nonnegative fields use sequential scales; correlation uses [-1,1]. Independent
instantaneous limits are common across all three representative frames and all
conditions. Representative frame keys 1, 879 and 1757 are the first, middle and
last numeric keys, not phase-selected samples or typicality claims. Quivers are
separate figures with a shared arrow scale and every fourth grid node shown.

## Reuse for another experiment

Add input paths and identifiers to `conditions` in the JSON configuration.
Provide documented coordinate_scale, velocity_scale and units. Supply an
orthogonal `transform` that maps the native basis to transverse x and axial y;
it is applied identically to coordinates and components. `origin` is expressed
in scaled native coordinates before transformation. Set mapping_verified only
after checking geometry. Set valid_flag_values only with a recorded
flag_definition_source. If physical velocity conversion needs calibration
divided by pulse separation, enter that verified combined velocity_scale and
record the evidence in notes. Unit scale defaults are diagnostic identity
operations and are not physical calibration.

The reader currently targets this HDF5 layout; adapt the reader for another
layout rather than silently guessing. New conditions are processed frame by
frame. Results, color limits, figure manifest and software versions are saved.

Validation includes independent pointwise means/RMS/covariance, nonnegative
variance, correlation bounds, the mean-speed inequality and grid shape checks.
No uncertainty, carrier-fluid stress or turbulence claims follow from these
uncalibrated statistics.
