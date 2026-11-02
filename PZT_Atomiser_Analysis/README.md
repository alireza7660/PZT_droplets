# PZT atomiser final analysis

The active configuration is `config/analysis.json`. Final deliverables are in
`Final_Results/Batch4_Cam1/`, with figures, data, diagnostics and methods
subdirectories. Earlier results in `provisional_calibrated/` retain the previous
5.05 micrometres/pixel calibration and remain archived unchanged. The previous
README is preserved as `README_5p05_archive.md`.

## Regenerate all final outputs

From the PZT_droplets workspace root:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.pipeline --config PZT_Atomiser_Analysis/config/analysis.json
```

The executable notebook is `notebooks/PZT_final_analysis.ipynb`; it calls the
same entry point. `notebooks/Spatial_maps.ipynb` is also updated to the final
workflow. To build and execute the notebook from the command line:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.build_final_notebook
```

## Calibration and Cartesian names

The confirmed PZT calibration is now **5.35 micrometres/pixel**, superseding both 5.05 and 7.25.
Pulse separation remains 25 microseconds. Raw us/vs pixel displacements are
converted once by 0.214 m/s before any statistics are calculated. Coordinates
use 0.00535 mm/pixel. No already calibrated result is multiplied by 0.214.
Original HDF5 files are read only.

Cartesian x/u is horizontal/transverse and positive right; y/v is vertical/
axial and positive upward. Existing array layout, component signs and the
exported (0,0) image-reference origin are retained. No array is transposed,
rotated or additionally reflected. The origin is not an inferred outlet.
The adopted export/component convention is stated once in the final methods.
No PZT calibration change is applied to any ZNMF source or result.

## Full-field signed contours

Mean and instantaneous axial-velocity contours, the station-location map and
flag-specific axial diagnostics show every finite signed measurement. Masks
apply only to genuinely missing or invalid data. Blue-white-red maps negative
values to blue, zero to white and positive values to red. Each mean component
uses symmetric limits from its own maximum absolute ensemble value, excluding
flag diagnostics. One linear normalization is used across the full range;
positive and negative values are not stretched separately. Contour levels are
uniformly spaced and signed fields include zero exactly. Near-zero values stay
near white. Comparable instantaneous fields share limits across the three
selected frames. No clipping, absolute values, filtering or gap filling is used.

Physical values are recalculated from the original displacements with the new
calibration. Component signs, spatial orientation, sample populations, station
indices, averaging and vector directions are preserved. Every
contour NPZ stores full values and a mask identifying only invalid data. Main
colour limits exclude flag-specific diagnostic extremes. The two in-plane
velocity magnitude definitions remain distinct supplementary quantities.

## Final deliverables and verification

- `figures/`: separate primary mean-v/mean-u contours, six other statistics,
  selected instantaneous fields/quivers, native profiles, station map and
  empirical PDFs/CDFs; vector PDF, SVG and 600-dpi PNG.
- `data/`: signed spatial statistics, instantaneous samples, exact station
  samples/coordinates, profiles, verification vectors and numerical figure NPZ.
- `diagnostics/`: retained counts, flag distributions and flag-specific
  population-sensitivity maps and data.
- `methods/`: exact configuration, methods, Markdown/LaTeX captions, figure/data
  manifest, software versions, validation, measured extrema and display limits.

Validation checks scaling against both the archived 5.05 calibration and the
saved 7.25 numerical reference. Velocities, RMS and coordinates scale by
5.35/5.05 or 5.35/7.25 respectively; variances and covariance by their squares;
correlations and counts remain unchanged. Station rows/columns and field keys
are retained. Numerical signed fields are checked independently of their
full-field validity masks. Methods and captions explain axes once;
figures use concise Cartesian mathematical labels without titles.

The superseded 7.25 statistics/configuration and its historical full-field
verification record are retained under
`Final_Results/Batch4_Cam1/diagnostics/calibration_reference_7p25/`.

Reusable entry points are `analysis/final_analysis.py`, the shared streaming
accumulator in `analysis/calibrated.py`, and `analysis/plotting.py`. The figure
manifest links every PDF/SVG/PNG to its signed numerical data source and caption.

## Fixed-x vertical axial profile

Full regeneration also finds the finite global maximum of mean_v and extracts
its native fixed-x column over the complete measured y range. Samples are
sorted by y together with their values, N and zero-based native indices.
Signed values and gaps are preserved; no adjacent-column averaging, smoothing,
interpolation or extrapolation is used. This is a fixed-x profile through the
global maximum rather than a tracked centreline.

To regenerate only the two added figures and the extracted profile from the
existing final statistics, without rerunning or changing other figures:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.vertical_profile
```

The figure stems are `mean_v_vertical_profile_location` and
`mean_v_vertical_profile`, with PDF, SVG and 600-dpi PNG in the final figures
directory. Extracted CSV/NPZ are in `data/mean_v_vertical_profile.*`.
Selection and verification metadata are in `methods/vertical_profile_*.json`.

## Plume-aligned measurement review ? interaction work paused

The current requested analysis is in
`Final_Results/Batch4_Cam1/Plume_Aligned/`. It fits one straight hydrodynamic
plume axis over a documented laboratory-height interval, checks sensitivity,
and rotates analysis copies of the native positions, velocity components,
covariance tensor and gravity using one orthonormal matrix. Contours use the
rotated native grid without resampling. Profile/section extraction retains
finite measured support without extrapolation. The original fixed-column
profile is preserved separately; the fitted-axis profile samples a different line.

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_plume_aligned_notebook
```

Configure the fitting interval and section stations in `config/plume_alignment.json`.
`transformation.json` documents angle, reference, matrix/inverse, gravity and
interval sensitivity; `data/extracted_centres.csv` records fitted centre residuals.
The origin is a measured-region analysis reference, not the atomiser surface.
The aligned axial coordinate y_p follows the plume and is not laboratory vertical.
Curvature, asymmetry and nonzero transverse velocity remain.

The active notebook shows only measured plume alignment at this stage. Previous
interaction notebook outputs are preserved separately in Plume_Aligned/validation;
the previous Combined_Flow_Prediction numerical results and figures remain unchanged.
`interaction_work_status=paused_for_plume_alignment_review` guards prediction
entry points; no interaction prediction or new interaction station is computed.
This coordinate transformation does not assert that physically remounting the
atomiser would reproduce the same flow. Original calibrated PZT data remain unchanged.

The completed aligned package uses x_p/y_p and a measurement-plane origin at the lowest supported laboratory height. Principal PDF/SVG/600-dpi PNG figures are under Plume_Aligned/figures/report; full native support is under diagnostics/full_support; before/after and count/fit diagnostics are under diagnostics/figures. Superseded generated packages are archived under Plume_Aligned/Archive. See its README and transformation.json for exact crop bounds and inverse mapping.
