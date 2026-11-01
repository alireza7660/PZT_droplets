# Fixed-line droplet mean velocities and stationary-air references

This package distinguishes a measured ensemble-mean profile from individual
monodisperse reference trajectories. Assumed diameters are [5, 10, 20, 50, 100]
micrometres. They are not measured or fitted sizes. Air at 20 C is recorded;
stationary air, pressure 101325 Pa and equal water temperature are assumptions.
Humidity is unavailable; dry-air properties are approximate and evaporation is omitted.

## Measurement and references

The native fixed line is x=14.54665 mm. It is not a tracked centreline.
Height h=y-y0 starts at y0=0.0856 mm (image reference), with measured
initial v=1.65501079 m/s and u=0.161064719 m/s. The location
is selected by finite paired means and N>=2; this
is sample support, not independent physical validation or an atomiser-surface location.
Every profile row, signed value, gap, native index and original y coordinate is retained.
Existing final calibration and sample selection are unchanged; SI velocities are not recalibrated.

Gravity, buoyancy and spherical finite-Re drag define the stationary-air references.
The same physical equations, properties and measured initial state used before this
revision are retained. Each upward branch stops at v=0 or the measured upper boundary.
`reference_endpoints.csv` separates actual turning heights from plotted endpoints.
Gravity-only actual rise is 139.65323 mm;
its plotted curve stops at the 27.2208 mm measurement
boundary, not at its turning point. The full and 0-2 mm figures contain identical
curve arrays and styles, with different horizontal viewing limits.

`reynolds_measured_velocity_assumed_air.csv` uses the magnitude of the measured
mean in-plane vector with stationary air and assumed diameters. This is not mean
instantaneous Re. It spans the measured domain. `reynolds_reference_trajectories.csv`
uses each model trajectory's instantaneous vector and ends at its own endpoint.
Gas density and viscosity enter both definitions; the missing third component is not invented.

`particle_relaxation_time.csv` gives tau_p,Stokes=rho_water*d^2/(18*mu_air).
Each trajectory CSV also gives tau_drag_secant=tau_p,Stokes/f(Re), which is an
instantaneous drag/slip coefficient timescale, not an exact finite-Re relaxation time.

The hypothetical ZNMF-based Stokes comparison is isolated in
`Jet_Interaction_Sensitivity/`. Its own downstream coordinate is preserved.
It is not simultaneous PZT gas data, a nozzle-exit Stokes number or evidence of entrainment.

Differences between measurements and stationary-air references do not by themselves
identify a coding error or a droplet size. Gas motion, size/initial-velocity distributions,
and lateral plume displacement may contribute; no dominant mechanism is established.
No unknown gas field or diameter distribution is fitted and no experimental model
validation is claimed.

## Exact report files

`report_manifest.json` lists the selected files. Use `methods_results.tex` for the
concise scientific section and `captions.tex` for standalone captions;
`report_figures.tex` supplies matching figure environments. These snippets require
amsmath and graphicx in the host report. Add references.bib to the host report's
bibliography. Figure paths are relative to this folder.

- `velocity_reference_full.pdf` (also 600-dpi PNG and SVG).
- `velocity_reference_near_origin.pdf` (also 600-dpi PNG and SVG).
- `reynolds_measured_velocity_assumed_air.pdf` (also 600-dpi PNG and SVG).
- `reynolds_reference_trajectories.pdf` (also 600-dpi PNG and SVG).
- `particle_relaxation_time.pdf` (also 600-dpi PNG and SVG).
- `diameter_sensitivity.pdf` (also 600-dpi PNG and SVG).
- `Jet_Interaction_Sensitivity/znmf_local_nozzle_stokes.pdf` (also 600-dpi PNG and SVG).

Every curve has CSV provenance in `figure_curve_manifest.json`; figures are drawn
directly from those CSV values. Per-diameter trajectory CSVs include original y,
relative height h, time, vector components, in-plane Re and defined timescales.
`diameter_sensitivity.csv` reports actual versus displayed endpoints; the full
endpoint table also includes gravity-only motion. `parameters_used.json` records
all assumptions, properties, source URLs, calibration, initial state and ZNMF provenance.

## Regeneration

From the PZT_droplets workspace root, run this single command using the project environment:

```powershell
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.update_theory_notebook
```

This rebuilds the theory package from final saved SI data, refreshes only the
theory section of the existing notebook and runs saved-output validation. Existing
experimental cells, figures and data are preserved. Configure assumptions in
`../../../config/theory_comparison.json`. Reusable physics is in
`../../../analysis/theory_comparison.py`, orchestration in `theory_package.py`,
and scientific text in `theory_report.py` in the same analysis directory.

## Review history, archive and validation

The local Git checkout had no commits later than reviewed commit 1522f54161c7.
The active theory source/configuration matched that checkout before editing.
`Archive_Before_Report_Revision_1522f54/` contains the superseded theory outputs
and source/configuration snapshots; it is excluded from report_manifest.json.
No final experimental output outside Theory_Comparison was archived.
Detailed signs, units, convergence, endpoints, calibration, figure/CSV agreement
and experimental-file hash checks are kept separately under `validation/`.
