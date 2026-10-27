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

The confirmed PZT calibration is now **7.25 micrometres/pixel**, replacing 5.05.
Pulse separation remains 25 microseconds. Raw us/vs pixel displacements are
converted once by 0.290 m/s before any statistics are calculated. No already
calibrated result is multiplied by 0.290. Original HDF5 files are read only.

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
values to blue, zero to white and positive values to red. Limits follow the
actual data range without symmetric endpoints. Signed contour levels include
zero exactly. Comparable instantaneous fields share limits across the three
selected frames. No clipping, absolute values, filtering or gap filling is used.

Calibration, component signs, Cartesian coordinates, sample populations,
statistics, profiles, distributions and vector arrows are unchanged. Every
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

Validation checks scaling against the archived 5.05 calibration: velocities,
RMS and coordinates scale by 7.25/5.05; variances and covariance by its square;
correlations and counts remain unchanged. Station rows/columns and field keys
are retained. Numerical signed fields are checked independently of their
full-field validity masks. Methods and captions explain axes once;
figures use concise Cartesian mathematical labels without titles.

Reusable entry points are `analysis/final_analysis.py`, the shared streaming
accumulator in `analysis/calibrated.py`, and `analysis/plotting.py`. The figure
manifest links every PDF/SVG/PNG to its signed numerical data source and caption.
