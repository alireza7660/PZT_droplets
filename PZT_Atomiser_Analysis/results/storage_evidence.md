# Storage-convention evidence audit

The supplied physical calibration is confirmed: 5.05 micrometres/pixel,
25 microseconds between PIV pulses, spray upwards, axial U positive upwards,
transverse V positive rightwards. The missing evidence is the original
processing/export definition, not the physical calibration.

## Workspace inspected

The only original file is `RAW_data/Batch4_Cam1.h5` (1,082,569,560 bytes).
There are no original processing scripts, project settings, software source,
HDF5 writer, images or calibration sidecars in this workspace, including its
hidden entries. `.venv` is the installed Python environment and
`PZT_Atomiser_Analysis` contains this session's generated analysis. In
particular, `config/analysis.json` is not the configuration used to produce
the HDF5 file.

`analysis/audit_storage.py` records the original-file inventory and visits all
HDF5 objects. Its output is `results/storage_evidence.json`. The earlier full
dataset inventory is `results/input_inventory.json`.

## Conclusions and exact evidence

| Question | Evidence in RAW_data/Batch4_Cam1.h5 | Conclusion |
| --- | --- | --- |
| Displacements or velocities? | `/PIV FIELDS/<key>/u`, `/v`, `/us`, `/vs` are float32 arrays without unit attributes. `/PARAMETERS/scale=1`, `/dt=1`, `/calib=0`, `/nondim=0`, `/Uc=1`, `/Lc=1`; no software definitions are stored. | Unresolved. These scalar names/values do not establish whether conversion happened upstream or what a calibration flag means. |
| Horizontal or vertical component? | The component names are u/v/us/vs with shape (319,80); coordinates are `/PARAMETERS/x` length 80 and `/PARAMETERS/y` length 319. | Grid array dimensions are established, but a component name or dimension is not proof of physical motion direction or sign. |
| Existing scaling, axis swaps or reflections? | Native x=31..5087 with step +64; native y=5104..16 with step -16. `/PARAMETERS/Xorigin=0` and `/Yorigin=0`. No transformation history, exporter name/version or coordinate units are stored. | Stored y decreases with row. This establishes array ordering only, not the image-to-physical transform or whether v was reflected consistently. No previous velocity scaling can be excluded. |
| Flag meanings? | `/PIV FIELDS/<key>/flag` is int32; full-record code counts are 0:1,097,567 and 1:43,741,073. There are no code definitions, enum attributes, validation settings or replacement-policy metadata. | Neither code can be labelled accepted, rejected, masked or interpolated from this evidence. |

There are no attributes on the root or any HDF5 object. The three root groups
are PARAMETERS, PIV FIELDS and STATISTICS. No textual dataset or userblock
provides software provenance. The complete numeric parameters are preserved
in the JSON audit.

An independently computed consistency check found that `/STATISTICS/U` and
`/STATISTICS/V` exactly match all-frame means of `us` and `vs`, respectively.
They do not match the corresponding u/v means: maximum absolute differences
are 0.1579631398 and 0.5165216043 native units. See
`results/Batch4_Cam1_metadata.json` and `analysis/pipeline.py` lines 102-108.
This establishes a numeric relationship, not whether us/vs mean smoothing,
replacement, scaling or another operation. No such semantics are inferred.

## What our code has done

`analysis/pipeline.py` lines 13-27 reads the selected components, applies the
configured orthogonal transform and velocity scale, and implements flag
selection only with a documented source. Current `config/analysis.json`
selects u/v, identity transform, null velocity/coordinate scales and no flag
selection. Null scales execute an identity factor. These are explicitly
diagnostic choices, not evidence about how the original file was produced.
Coordinate handling is in `analysis/pipeline.py` lines 35-40. No additional
calibration or sign reversal has been applied by this analysis.

## Precisely what to supply

Provide the **original processing/export project folder that produced
Batch4_Cam1.h5**, including:

1. The script/module/function that creates the HDF5 datasets named
   `PIV FIELDS`, `PARAMETERS`, `STATISTICS`, `u`, `v`, `us`, `vs`, and `flag`.
   Include called functions that calculate or modify those arrays before
   writing. The writer may use h5py File/create_dataset, PyTables, MATLAB
   h5create/h5write, or another HDF5 API; the dataset names are the most
   specific search terms available from this file.
2. The actual saved settings/project or batch-launch script used for Batch4
   Cam1, including any spatial/time conversion, component selection,
   transpose/swap/reflection and export options. Generic software defaults
   alone cannot establish which options were used in this run.
3. The validator/flag enum definition and the replacement/interpolation code
   assigning or modifying flag, u/v and us/vs. This may be in the same source
   tree as the writer.
4. One original image or documented image-coordinate convention linking the
   processing crop/rotation to the export. An outlet reference (c0,r0) is
   needed only if an outlet-centred origin is desired; a documented image
   reference is otherwise sufficient.

If the writer is proprietary, provide the software name and exact version,
the saved Batch4 Cam1 project/export-settings file, and that version's HDF5
format and flag-definition documentation. The exact original filename or
folder path cannot be recovered from this HDF5: none is embedded. Inventing
one would make the request less precise. Copy the above original project
folder into this workspace, for example as `Processing_source/` (a suggested
destination, not a claimed existing folder).

## Conversion decision

If source proves components are pixel displacements per pair, apply the
single factor (5.05e-6)/(25e-6)=0.202 m/s per pixel. If source proves they are
already in m/s, use factor 1. If another unit convention is established,
derive its conversion from that convention. Apply the verified component
mapping to coordinates and velocities together. Downward-positive image
vertical displacement needs a minus sign for upward U only if the exporter
has not already applied it. Pixel-coordinate conversion is 0.00505 mm/pixel
only if the stored coordinates are confirmed pixels.

Until the original exporter/settings are available, none of these branches
can be chosen defensibly. No new figure package was generated for this audit.
