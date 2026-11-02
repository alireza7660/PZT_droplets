# Complete plume-aligned measurement package

Interaction predictions remain paused. No crossing location is selected or updated.

## Reproduction
From the local PZT_droplets project root:
```powershell
.\.venv\Scripts\python.exe -m PZT_Atomiser_Analysis.analysis.update_plume_aligned_notebook
```
This executes the alignment notebook section, rebuilds all outputs and saves the executed notebooks. Existing generated outputs are archived under Plume_Aligned/Archive/<timestamp> before replacement. Configuration: config/plume_alignment.json. Source: analysis/plume_alignment.py and analysis/plume_alignment_complete.py.

## Measurements and axis estimate
Already calibrated authoritative laboratory fields are used unchanged: 5.35 micrometres/pixel, 25 microseconds, 0.214 m/s per pixel/pair. No calibration is applied here. Original finite-pair sample counts and masks are preserved. Original images linked to this PIV sequence were unavailable; previously inspected nearby EVOS microscope images and assembly illustration do not establish acquisition tilt. The archived image_inventory.json documents that inspection.
One subgrid quadratic peak of laboratory mean_v is extracted per resolved native y row; equal-weight least squares fits x_c=a+b*y. Fit angle: 5.259017932 degrees right of laboratory vertical. Actual fit interval: 5.0504 to 21.9992 mm, 199 centres. RMS horizontal residual 0.0326113 mm; maximum absolute residual 0.0790231 mm. Interval sensitivity: 4.864496 to 5.431601 degrees. This is a hydrodynamic straight-axis estimate, not independently established mounting tilt. Residual curvature and asymmetry are not attributed entirely to mounting. Subgrid localisation does not imply corresponding instrument accuracy. Secondary half-maximum centroid fitting and interval results are recorded in transformation.json.

## One rigid transformation and measurement reference
Origin in laboratory coordinates: (14.510910931, 0.085600000) mm, the fitted-axis intersection with the lowest supported laboratory-y measurement plane. It is the measurement reference plane, not the physical atomiser surface. x_p is transverse, y_p increases along the plume; y_p is not laboratory vertical.
Position p_p=Q(p_lab-p_ref); velocity w_p=Q*w_lab; covariance C_p=Q*C_lab*Q^T. Inverse p_lab=Q^T*p_p+p_ref and w_lab=Q^T*w_p. Q and its inverse, original lab coordinates, calibrated lab source paths and transformed gravity [0.898861356725992, -9.765369039820497] m/s^2 are exported. RMS values are square roots of rotated tensor diagonal entries, never rotated as a vector. No rows are shifted independently, mirrored, averaged or symmetrised; nonzero transverse velocity and actual curvature remain. This representation does not establish what physical remounting would produce.

## Native support and report window
Full-support native-grid maps: figures/full_support/. Principal standalone report figures: figures/report/. Before/after, fit, residuals, sample counts and original-column comparisons: diagnostics/figures/.
Report x_p bounds [-5.5, 5.5] mm; y_p bounds [0.5062519851358152, 26.829617988796766] mm. Width covers all extracted half-maximum lobes plus the configured margin, rounded upward; given this width, the longest axial interval with all four corners in measured laboratory support is retained. Retained axial span 96.296% of the fitted-axis measured span and 88.974% of full rotated bounding-box vertical extent; excluded 3.262170 mm. Native nodes inside the report rectangle: 38.70%. This node fraction includes deliberate lateral crop, not data rejection. All native nodes remain exported and shown in full-support maps. Report windows clip contour polygons; no window stretching, missing-corner filling or extrapolation occurs. Internal missing-data masks are preserved (corner_mask=False). Spatial aspect ratio is equal and transverse display bounds symmetric. Signed maps use symmetric linear blue-white-red normalization with uniformly spaced levels; matching full/report maps share limits, instantaneous components share limits across all selected frames. RMS starts at zero; correlation spans [-1,1].

## Profiles and instantaneous measurements
The fitted-axis profile samples the fixed straight x_p=0 line at native laboratory y planes, using four-finite-corner bilinear interpolation without extrapolation. It differs from the original fixed laboratory-column line, whose original exports are copied byte-for-byte to reference/. Transverse cuts are fixed y_p stations [2, 5, 10, 15, 20, 25] mm; source lab positions, signed components and finite support are saved. Full section values are exported even where the report display clips lateral tails. No smoothing or adjacent-column averaging is used.
Selected saved instantaneous frames are discovered from original data/frame_*.npz (currently 1, 879, 1757). Actual already-calibrated u/v arrays are rotated, not reconstructed from means. Both component maps and separate quivers are exported in full support and report form. Quivers subsample displayed arrows only and share a physical scale.

## Verification and retained asymmetry
Inverse position/vector transforms, vector-magnitude preservation, covariance inverse/trace/eigenvalues, unchanged native N and masks, outside-support sampling, and direct profile-plot/export agreement are checked. Original laboratory numerical data and figures are protected by SHA256 comparisons. Instantaneous inverse/magnitude/mask checks are recorded separately. No interaction integration is run.
Station metrics in data/transverse_profile_metrics.csv retain left/right half-width differences; genuine asymmetry and nonzero transverse velocity remain. These descriptive differences are not a formal statistical significance test.

## Outputs
38 separate figures, each PDF, SVG and 600-dpi PNG. report_manifest.json enumerates exact paths and captions; numerical native grids, original lab coordinates, profiles, full tensors, instantaneous arrays, support and transformations are in data/ and transformation.json.
