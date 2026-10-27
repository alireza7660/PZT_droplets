# PZT atomiser methods

Coordinates are Cartesian: $x$ is horizontal/rightward and $y$ vertical/upward; $u$ is transverse and $v$ axial. The in-plane vector is $\mathbf{w}=(u,v)$; primes denote fluctuations about the respective means.

Calibration is 7.25 micrometres/pixel (0.00725 mm/pixel), replacing 5.05 for
this PZT experiment only. Pulse separation is 25 microseconds. The adopted
export convention is us/vs as horizontal/vertical pixel displacements per
pair, with retained signs: u=0.290 us and v=0.290 vs (m/s). The original HDF5
is opened read-only, converted once at the instantaneous-sample stage, then
all statistics are calculated. Naming is changed without transposing,
rotating or reflecting arrays. The exported (0,0) image reference remains
the spatial origin; it is not an established atomiser outlet centre.

The finite-pair population requires both components finite. No numeric flag
code is rejected or assigned a physical validity meaning. RMS, variance and
covariance use denominator N. Undefined correlation and fewer-than-two-pair
statistics are masked. No derivative, jet-width or similarity analysis is
added. Upstream replacement is not inferred. No additional spatial smoothing,
gap filling or extrapolation is applied.

All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints. This presentation applies to mean and instantaneous axial
contours, the station-location map and flag-specific mean-axial diagnostics.
Profiles, PDFs/CDFs, quivers and numerical arrays also retain signed values.
All component velocities, covariance and correlation use signed colour scales.
Main limits are calculated from the actual ensemble and selected instantaneous
frames; flag-specific extrema do not set ensemble colour limits. Diverging
normalization centres zero consistently without imposing symmetric endpoints.
Signed contour boundaries include zero exactly. Only cells touching genuinely
missing or invalid nodes are omitted (corner_mask=False); finite negative
measurements do not produce artificial cutouts.

Magnitude of the mean vector is sqrt(mean(u)^2+mean(v)^2). Mean instantaneous
in-plane velocity magnitude is mean(sqrt(u^2+v^2)); the definitions are distinct.
Mean component contours are primary; magnitude quantities are supplementary.
Native station indices and three field keys are preserved. Profiles use the
same native rows without interpolation; PDFs are separate native-node histogram
densities with unit integral and CDFs are exact empirical staircases. No frame
independence, turbulence or absolute-uncertainty claim is made.

Grid spacings are 0.464 mm horizontally and 0.116 mm vertically; these describe
vector spacing, not interrogation-window size or effective resolution.

Regenerate from the PZT_droplets workspace root:
    .\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.pipeline --config PZT_Atomiser_Analysis/config/analysis.json
