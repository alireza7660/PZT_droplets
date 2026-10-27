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

## Axial contours: display policy only

Mean and instantaneous axial-velocity contour maps display only values >= 0,
with a sequential colour scale starting at zero. Negative nodes are masked,
not replaced with zero or absolute values. Cells touching a masked node are
omitted without gap filling. Every contour NPZ contains the full signed values
and a separate display_mask. This also applies to the station-location map and
flag-specific mean-axial contour diagnostics. Captions state the display choice.

Statistics, numerical arrays/CSV, PDFs/CDFs, profiles and vector arrows retain
all signed component values. Transverse velocity, covariance and correlation
retain signed colour scales. Mean in-plane velocity-vector magnitude and mean
instantaneous in-plane velocity magnitude remain distinct supplementary results.

The finite-pair population requires both us and vs finite, without flag-code
rejection. Flag counts and code-specific mean/RMS diagnostics remain separate
and do not set main-ensemble colour limits. Main mean-v colour limits start at
zero and reach the actual nonnegative maximum. Other signed field limits use
actual finite ensemble values and a zero-centred diverging normalization.
Selected instantaneous fields share colour limits. No filtering, smoothing,
extrapolation, gap filling or new jet-fit/width analysis is introduced.

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
nonnegative axial display masks. Methods and captions explain axes once;
figures use concise Cartesian mathematical labels without titles.

Reusable entry points are `analysis/final_analysis.py`, the shared streaming
accumulator in `analysis/calibrated.py`, and `analysis/plotting.py`. The figure
manifest links every PDF/SVG/PNG to its signed numerical data source and caption.
