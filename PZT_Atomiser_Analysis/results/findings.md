# Findings and measurement limits

The HDF5 contains 1757 fields, each with 319 rows and 80 columns. Every native
u/v pair is finite. Population moments were independently checked at three
grid locations; the code also checks correlation bounds, nonnegative variance
and mean instantaneous speed >= magnitude of the mean vector.

For the paired-finite native u/v population:

| Quantity | Spatial minimum | Spatial maximum |
| --- | ---: | ---: |
| Mean native u | -0.117344 | 1.024717 |
| Mean native v | -0.119483 | 7.995459 |
| Magnitude of mean vector | 0.000429 | 8.035626 |
| Mean instantaneous speed | 0.089260 | 8.966638 |
| Native u RMS | 0.067266 | 3.608221 |
| Native v RMS | 0.065353 | 4.266189 |
| Signed cross-covariance | -0.638639 | 2.496695 |
| Correlation (dimensionless) | -0.957429 | 0.959849 |

Component and speed units remain unverified native units; covariance units
are their square. These extrema are grid-node extrema over full-record
pointwise moments, not pooled samples or interpolated extrema. Mean native v
is largest at native x=2719, native y=16, on the measurement boundary. Its
maximum therefore does not establish a physical outlet velocity or the global
spray maximum outside the field of view.

Flag code 0 occurs in 2.4478% of all stored node-frame entries, with local
fractions up to 51.0529%. These are code frequencies; acceptance, interpolation,
droplet concentration and intermittency cannot be assigned to them. All
node-frame entries remain finite, so an all-finite fraction map is uniformly
one and should not be mistaken for a successful quality-validation result.

Stored statistics U/V equal means of us/vs exactly. They differ from native
u/v means by maximum absolute differences of 0.157963 and 0.516522,
respectively. Mixing stored moments with recomputed u/v speed or covariance
would therefore mix populations/processing definitions. All current figure
statistics are consistently recomputed from u/v. No interpretation of the
meaning of us/vs or the upstream processing follows from this comparison.

The user confirmed 5.05 micrometres/pixel, pulse separation 25 microseconds,
upward-positive axial U and rightward-positive transverse V. The factor
0.202 m/s per pixel is conditional on pixel-displacement storage. The native
y array already decreases with row index, so blindly reversing the vertical
component could reverse an existing export transformation. No export code,
software identification, original image, flag definition or outlet pixel
reference is available. Physical calibration, component mapping and quality
selection cannot yet be completed defensibly. The pipeline records the
confirmed settings and keeps these diagnostics explicitly uncalibrated.
