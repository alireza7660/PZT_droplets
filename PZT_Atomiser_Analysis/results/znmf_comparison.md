# Targeted PZT / established ZNMF workflow comparison

Physical PZT calibration remains confirmed: 5.05 micrometres/pixel and 25
microseconds. Spray direction is upward; physical axial U is upward positive
and transverse V rightward positive. None of those facts establishes how the
exporter stored its components.

## Files directly compared, read only

- PZT: `PZT_droplets/RAW_data/Batch4_Cam1.h5`.
- ZNMF: `ZNMF_Results/Raw data/First_location/Batch4_Cam1.h5`.
- ZNMF: `ZNMF_Results/Raw data/Second_location/Batch4_second_Cam1.h5`.
- ZNMF: `ZNMF_Results/Raw data/Third_location/Batch4_third_Cam1.h5`.

These paths share the parent `Inhaler Project`. The comparison implementation
is `analysis/compare_znmf.py`; full metadata and comparisons are saved in
`results/znmf_comparison.json`. No ZNMF file was modified, no processing
workflow was executed, and no figures were generated.

## Direct HDF5 comparison

| Feature | PZT | ZNMF, all three raw FOVs |
| --- | --- | --- |
| Root groups | PARAMETERS, PIV FIELDS, STATISTICS | Identical |
| Field names/types | flag int32; u,v,us,vs,omega float32 | Identical |
| Statistics names/types | NoSamples int64; U,V,uu,vv,uv,Omega,omega_sq float64 | Identical |
| Parameter names | Same full set | Identical |
| Calibration-like values | calib=0, nondim=0, scale=1, dt=1, Uc=1, Lc=1 | Identical |
| Origin parameters | Xorigin=0, Yorigin=0 | Identical |
| Other scalar parameters | ny=319, stats=1, viscosity=1 | Identical |
| nx / fields | 80 / 1757 | 79 / 4000, 3850, 3780 |
| x coordinates | 31 to 5087, step +64 | 63 to 5055, step +64 |
| y coordinates | 5104 to 16, step -16 | Identical |
| Attributes, including units/software/version | None anywhere | None anywhere |

Only nx and NumberOfFields differ among scalar PARAMETERS. The schema is
strongly compatible with a shared export format. It is not proof that the
same software version, calibration mode or component transformation was used.
The x-grid offset/extent also differs, so identical sampling/export settings
cannot be asserted from these arrays.

## What the established ZNMF workflow actually documents

1. `ZNMF_Results/plot_mean_velocity.py:8-14` describes its input mean U as
   horizontal pixel displacement, right positive; it applies a negative
   conversion for the ZNMF leftward downstream direction and retains its y
   convention. That is documentation in an analysis reader, not the original
   HDF5 writer or its version/settings.
2. `ZNMF_Results/Final_Result/analysis/pipeline.py:125-138` calculates the
   velocity factor from calibration/timing and loads `us` and `vs`, with
   minus and plus signs respectively. Lines 160-161 construct the ZNMF
   coordinates. These ZNMF-specific signs, origin, nozzle scale and FOV
   registration are not transferred to PZT.
3. `ZNMF_Results/Final_Result/report/methods.md:65-76` explicitly states
   `U=-0.202 us`, `V=+0.202 vs` for ZNMF, finite paired acceptance, and that
   flag meanings are not established. The matching array/statistics
   relationship is also observed in PZT: stored U/V match means of us/vs.
   This supports a common numerical relationship, not documented semantics
   for u/v versus us/vs or pixel-versus-SI storage.
4. `ZNMF_Results/numerical_validation_audit/numerical_audit_report.md:13`
   explicitly lists exporter flag meanings and alternate u/v versus us/vs
   processing as unresolved. Lines 58-64 distinguish uncalibrated-looking
   metadata and the existing pixel-displacement interpretation from
   independent exporter certification.
5. `ZNMF_Results/Final_Result/analysis/pipeline.py:118,481` records unknown
   flag semantics and accepts finite us/vs pairs. The earlier
   `ZNMF_Results/pdf_velocity_analysis.py:195,222,255` selects flag 1 but
   explicitly acknowledges that its physical meaning is undocumented.
   Thus flag 1 is not a convention validated by the established workflow.

The inspected ZNMF reports/scripts do not identify the original PIV software
name/version or saved export settings. Camera_setup.png was inspected: it is
a FOV-registration schematic, not a processing-software settings screenshot.
The sibling PIV folder contains a reference PDF, not an exporter installation
or settings project. The PZT_droplets_GitHub folder contains copies of this
session's analysis files, not original processing source.

## Supported reuse and remaining decision

Reuse the reader schema, numeric field ordering, actual coordinate arrays,
and the verified distinction between stored-statistics components and the
alternate pair. The scalar-parameter and coordinate comparison supports a
pixel-export hypothesis but cannot exclude different upstream conversions.
Do not transfer the ZNMF sign reversal, origin, nozzle diameter, registrations
or flag-1 selection. Do not substitute us/vs for u/v until the processing
relationship, including filtering or replacement, is explained.

The factor 0.202 is exact for the supplied calibration and timing **if** a
stored component is pixel displacement per pair. Prior SI velocity requires
factor 1 instead. Matching formats alone does not choose between those
branches. Physical regeneration remains pending storage and component
definitions; no duplicate uncalibrated package is needed.

## Targeted evidence request: source code is not required

Any screenshot, saved-setting excerpt or software-documentation entry that
establishes the following is acceptable:

1. **Software/version and run identity:** Help/About showing application and
   exact version, plus a project/result/export-settings view identifying
   Batch4 Cam1 for PZT and the ZNMF reference. If the same installation was
   used, explicitly confirm that both exports used that version. A current
   About screen alone does not establish the historical export version.
2. **Stored quantity and transformations:** the export/result settings for
   the HDF5 arrays showing displacement [pixel/pair] versus velocity [unit],
   calibrated/raw export mode, and any rotation, component swap, x/y flips
   or sign reversal. The view or corresponding HDF5-format manual entry
   must identify u/v and us/vs, their horizontal/vertical positive directions,
   and whether us/vs have been smoothed or replaced. Also establish whether
   PARAMETERS/x,y are pixels or physical coordinates. If no single settings
   page shows this, a format-manual definition plus run-settings screenshot
   is sufficient. Do not infer signs from the visible spray direction alone.
3. **Flag legend/replacement policy:** the version-specific manual/help entry
   or validation settings explicitly defining stored codes 0 and 1 and how
   rejected or replaced vectors are represented. It must say whether codes
   describe pre-replacement or post-replacement vectors and which arrays
   contain replacements. A coloured display without the numeric-code legend
   is insufficient.

An outlet pixel reference is optional if a clearly labelled image-reference
origin is used; no outlet-centred origin will be invented. Once these entries
resolve the definitions, apply calibration exactly once, map components and
coordinates consistently and regenerate only the physical outputs.
