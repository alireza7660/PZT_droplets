# Jet/PZT comparison

This folder is independent of both original projects. No original results were modified. See caption.txt for the figure caption and verification.json for exact input hashes, retained ranges, counts and numerical checks.

## Final inputs

- ZNMF: `C:\Users\LabAdmin\OneDrive - Monash University\Desktop\Inhaler Project\ZNMF_Results\Final_Result\tables\FOV_1_2_3_centreline_and_width.csv`
- PZT: `C:\Users\LabAdmin\OneDrive - Monash University\Desktop\Inhaler Project\PZT_droplets\PZT_Atomiser_Analysis\Final_Results\Batch4_Cam1\data\mean_v_vertical_profile.csv`

ZNMF provenance: Final_Result/analysis/figures.py plots accepted tracked centreline Uc; analysis/saved_results.py checks the centreline CSV against the saved renderer arrays. The selected full-domain measured maximum is 2.680478380412314 m/s at x/D0 = 25.2096, outside this comparison. It is not marked. Only status=measured, finite Uc values within the interval are retained.

PZT provenance: analysis/vertical_profile.py writes the final CSV and NPZ used for mean_v_vertical_profile. Final_Results/Batch4_Cam1/methods/vertical_profile_selection.json confirms x = 14.54665 mm. configuration_used.json confirms the final 5.35 micrometres/pixel, 25 microseconds and 0.214 m/s per pixel displacement calibration. Saved mean_v is already calibrated and is used unchanged. The experiment is the PZT atomiser at 190 kHz and 9 Vpp.

## Coordinates and outputs

ZNMF x_over_D0 is multiplied once by D0 = 2 mm. The reusable script determines units from explicit column names: x_over_D0 gets factor 2; x_mm or s_mm gets factor 1. Ambiguous headers are rejected. PZT s_mm = y_mm, with no coordinate rescaling or velocity recalibration. Each curve is cropped by its own measured coordinates to 0 <= s_mm <= 27. Native order, values and origins are preserved; no zero or endpoint samples are added. These are respective analysis coordinates from separate experiments, without physical outlet registration.

Retained ZNMF: 75 points, 2.2624-26.8256 mm, maximum 2.070241610919579 m/s. Retained PZT: 315 points, 0.0856-26.964 mm, maximum 1.6865481421305646 m/s. Common axes: s = 0-27 mm; velocity = 0-2.4 m/s. CSV files include original source columns and zero-based source data-row indices for audit. No smoothing, fitting, interpolation or extrapolation is performed.

PDF and SVG are vector exports; PNG is 600 dpi (3840 x 2640 pixels). The script verifies line arrays and exact CSV numerical round trips against source selections, the full-domain ZNMF maximum, fixed PZT x, and unchanged source hashes.

## Regenerate

Requires Python 3, numpy and matplotlib. From this folder in PowerShell:

```powershell
python .\plot_comparison.py --znmf '..\..\ZNMF_Results\Final_Result\tables\FOV_1_2_3_centreline_and_width.csv' --pzt '..\PZT_Atomiser_Analysis\Final_Results\Batch4_Cam1\data\mean_v_vertical_profile.csv'
```

Both source paths can be replaced with absolute paths. Optional `--output-dir` selects a different destination.

## Export review

The exported 600 dpi PNG was visually inspected: labels and legend are readable, curves are unobscured, borders and inward major/minor ticks are present, with no title or maximum annotation. PNG dimensions and resolution metadata were checked. The PZT source CSV also matches the final figure NPZ and the indexed final statistics field exactly. A separate PDF raster review was unavailable because no PDF renderer was installed and the package index could not supply pypdfium2; all three formats were exported from the same verified Matplotlib figure. The PDF skill's operation marker could not run because Node.js was unavailable.
