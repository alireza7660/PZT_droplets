# Provisional calibrated results — not final validated experimental results

Spatial calibration 5.05 micrometres/pixel and pulse separation 25 microseconds
are confirmed by the user. Storage units, exported component definitions,
vertical sign and flag semantics are working assumptions, not independently
verified facts. The source arrays are us/vs because their means reproduce
exporter STATISTICS/U,V. Their upstream smoothing/replacement is undocumented.
No additional filtering, spatial smoothing, gap filling or extrapolation is
applied here.

For the current configuration, V=(0.202) us + (0) vs and U=(0) us + (0.202) vs in m/s. The factor
is applied exactly once to instantaneous arrays, and all moments are computed
from those samples. The stored vs sign is retained provisionally; upward spray
direction and decreasing exported y do not establish its sign. No ZNMF sign
reversal is inherited. Nearby original images lack an identifiable mapping to
these export vectors and do not independently resolve signs.

Coordinate convention: x is axial/upward, y is transverse/rightward; U and V
are the corresponding velocity components. Lowercase u' and v' denote their
fluctuations about the means. The instantaneous in-plane velocity vector is
bold u=(U,V). Its magnitude is sqrt(U squared + V squared). Magnitude of the
mean in-plane velocity vector is sqrt(mean(U) squared + mean(V) squared);
mean instantaneous in-plane velocity magnitude is mean(sqrt(U squared + V
squared)). These distinct quantities are supplementary; mean_U and mean_V
are the primary velocity contours in the notebook and report selection.

Horizontal plotted y is 0.00505 times exported PARAMETERS/x;
vertical plotted axial x is 0.00505 times exported PARAMETERS/y. Origin is the
exported (0,0) image reference, not an established outlet. Both coordinate
arrays retain alignment with native rows/columns. The scaled grid spacings
are 0.3232 mm horizontally and 0.0808 mm vertically; these are vector spacing,
not interrogation size or effective resolution.

The finite-pair population consists of samples with both us and vs finite.
No flag is selected or rejected. All moments, including velocity magnitude, covariance
and component RMS, use this paired population and population denominator N.
At least two pairs are required for reported moments. Correlation is masked
at zero variance. N and finite-pair fraction are not vector-validation or
concentration metrics. Code-specific counts/means/RMS are selection-sensitivity
diagnostics; no physical acceptance or rejection meanings are assigned.

Profiles follow selected native rows without shifting or interpolation.
Stations are configured exact nodes, not pooled spatial ROIs or a physical
centreline. Locations are shown on sampling_stations and in stations.csv.
Histogram PDFs use common bin edges per quantity across stations, full sample
support and unit integral. CDFs are exact empirical staircases. Normalized
fluctuations use local N-denominator RMS; zero variance is explicitly omitted.
No KDE or negative-velocity magnitude extrapolation is used. No independence-based
confidence intervals, turbulence interpretation or periodic/stochastic
decomposition are claimed. Timing beyond the pulse separation and tracked
phase remain unverified. First/middle/last numeric field keys are illustrative,
not phase-selected or certified typical fields.

Nonfinite samples, low-count statistics and zero-variance correlations are
masked. Contours use native quadrilaterals and corner_mask=False; intrinsic
contour-boundary interpolation does not increase measurement resolution.
Common full-range limits include flag populations when comparing component
mean/RMS, and all conditions. Instantaneous scales are common across selected
frames/conditions. Every figure has PDF, SVG, 600-dpi PNG and numerical NPZ;
spatial tables, station samples and verification vectors preserve underlying
values. No plot titles are used.

To revise definitions, edit config/analysis.json and rerun the notebook or CLI.
components chooses the source pair. component_mapping maps it into [V,U];
coordinate_mapping maps exported coordinates into [transverse,axial]. transform
is a common orthogonal basis change applied to BOTH coordinates and velocities.
A later basis reflection therefore also changes arrows and recomputed signed
covariance consistently. If evidence identifies an exporter-only component
sign mismatch with already-reflected coordinates, correct component_mapping
and document why coordinate_mapping differs. Do not reflect an already-correct
coordinate twice. All statistics, samples, profiles and captions regenerate.
No derivatives are computed in this run; any future derivatives must use this
same transformed Cartesian basis. General rotations turn row profiles into
oblique cuts: the pipeline rejects such row-profile geometry rather than
mislabeling it as constant axial position.

## Configuration used

- us/vs are horizontal/vertical pixel displacements per image pair; not independently verified. Their means reproduce STATISTICS/U,V.
- Exported coordinates are pixels, scaled once to mm; outlet origin is not established.
- Stored vs sign is retained provisionally; spray direction and descending y do not verify its sign.
- Use the finite-pair population without selecting flag codes. Flag-specific populations are diagnostics, not validated vectors.
- us/vs may include upstream smoothing or replacement; this analysis adds none. Relationship to u/v remains unverified.

Effective coordinate matrix: [[1.0, 0.0], [0.0, 1.0]]; effective component matrix: [[1.0, 0.0], [0.0, 1.0]].
