# Prospective horizontal side-jet placements

Successfully integrated: 27 representative droplets,
across 9 PROPOSED arrangements. This is not measured combined geometry.
Jet centreline heights: [6.933599999999999, 13.696, 20.544] mm, selected at nearest native PZT rows
to 25%, 50% and 75% of the measured height range. Nozzle-to-plume distances:
[14.220799999999995, 31.996799999999986, 50.41919999999998] mm, selected at 20%, 50% and 80% of positive measured gas
centreline station indices. The measured local gas centreline, not its nominal
nozzle axis, is aligned at the PZT fixed-line anchor. Three adjacent native columns
41, 42, 43 initialise both measured signed components at each selected height.
No travel from an unmeasured atomiser surface is simulated.

All 27 cases reach the downstream measured-gas boundary before the 30 ms limit.
Endpoints are support limits, not droplet stopping or capture. Their displacements
occur over different elapsed times and must not be ranked as collection performance.
No collection target is specified and no optimum or capture efficiency is claimed.

Files: placement_comparison.csv (nine arrangements), case_comparison.csv (27 cases,
initial vectors, displacement and termination), H*_L*_C*.csv (27 full histories).
For each H1..H3 / L1..L3 pair, PDF and 600-dpi PNG figures are:
*_trajectories, *_velocity, *_displacement, *_particle_Re, *_supporting_Stokes.
Displacement components have separate ordinates; velocities retain signed values.
CSV histories also contain both gas components, vector slip, Re and defined response
timescales. Supporting Stokes uses local D0 / gas-vector magnitude, not centreline
Uc or a nozzle-exit timescale. See methods_results.tex for the prescribed-field
model approximation, source_manifest.json for sources, manifest.json for exact
selection/configuration and convergence.json for all half-step comparisons.

Regenerate with the single command in the parent README. To refresh notebook
display from existing computed CSV/figures without new integration, use
update_combined_flow_notebook --reuse-saved-studies. This is a display refresh only.
