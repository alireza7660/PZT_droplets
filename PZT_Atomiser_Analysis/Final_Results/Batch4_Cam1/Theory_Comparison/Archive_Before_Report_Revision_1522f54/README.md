# Droplet theory comparison: conditional references

Final PZT data: 190 kHz, 9 Vpp; 5.35 micrometres/pixel, 25 microseconds, 0.214 m/s per pixel displacement. Saved SI fields are read without applying calibration again. x/u are horizontal/rightward, y/v upward. All signed experimental velocities remain unchanged.

## Available evidence and missing inputs

`../data/statistics.npz` provides ensemble mean u/v, RMS, covariance and N; `../data/frame_1.npz`, frame_879 and frame_1757 provide representative instantaneous vector fields. RAW_data/Batch4_Cam1.h5 contains correlation-export fields us/vs and flags, not droplet identities or trajectories. analysis/calibrated.py accumulates finite paired samples; analysis/final_analysis.py adopts the final convention. No validated sizing result or simultaneous gas field was found in the active project or RAW_data. Image-spot width, the water-flow-rate image and nozzle dimensions are not diameter evidence. The exported (0,0) reference is not a located atomiser surface. The adopted component/displacement interpretation remains documented in the experimental configuration; this comparison does not independently verify it.

The fixed-x profile is column 42 through the global mean-v maximum, not a tracked centreline. A declining value there may partly reflect lateral plume displacement. Neither a fixed-position profile nor a height-wise maximum of ensemble means is an individual droplet trajectory. No height-wise maximum is substituted here.

## Properties and model

Air temperature is recorded as 20 C. Pressure 101325 Pa, water also at 20 C, and constant spherical diameter are assumptions. Humidity is unavailable; dry-air properties are an approximation, not a recorded zero-humidity condition. Evaporation is neglected. rho_air=1.20479 kg/m3; rho_water=998.206746 kg/m3; mu_air=1.81332212e-05 Pa s; g=9.80665 m/s2. Water density uses NIST GLP10's Tanaka polynomial; air uses the approximate NIST dry-air reference; Sutherland viscosity uses constants reported in NASA's wake-validation paper. No liquid viscosity enters Re or relaxation time.

Vector equations use dq/dt=w_d and dw_d/dt=(0,-g(1-rho_air/rho_water))+(w_air-w_d) f(Re)/tau_p, where Re=rho_air*d*norm(w_air-w_d)/mu_air, tau_p=rho_water*d^2/(18*mu_air). f=1 for Re<=0.01 (creeping-flow approximation); otherwise f=1+0.15 Re^0.687, equivalent to C_D=24*f/Re. This implementation enforces Re<800 following the cited COMSOL range. Using slip*f/tau avoids division by Re at zero slip. The rigid-sphere approximation omits internal circulation, deformation, evaporation, added mass, history forces, collisions and two-way coupling. These may limit quantitative applicability.

Both in-plane components enter drag and Re; the unmeasured out-of-plane component is omitted, so Re is an in-plane-slip estimate. Re of a mean vector is not the ensemble mean of instantaneous Re. No instantaneous sample population is relabelled as a measured diameter-resolved quantity. Assumed diameters [5, 10, 20, 50, 100] micrometres are sensitivity cases, not an inferred size distribution.

The initial state is the lowest finite paired mean-u/v sample with N>=2, y0=0.0856 mm, v0=1.65501079 m/s at fixed x=14.54665 mm. This is a sample-availability criterion, not independent physical validation. It is not an outlet launch condition. All PZT comparison heights use h=y-y0; absolute image-reference y values are also retained in CSVs. RK4 integrates individual monodisperse reference trajectories in assumed stationary air, including horizontal motion, and stops at the first v=0 turning point or measured upper boundary. Each reference starts at the measured x but may drift horizontally: it is not a prediction along the fixed-x measurement column. Gravity-only ballistic motion excludes both drag and buoyancy. No descending branch is joined to the ascending curve; no terminal settling prediction is compared to upward plume velocity.

## Response question and Stokes number

The particle relaxation time is defined as the Stokes-limit baseline tau_p=rho_water*d^2/(18*mu_air). It is reported for every assumed diameter, even when finite-Re drag is needed; it is not represented as an exact finite-Re relaxation time. Along trajectories tau_drag_secant=tau_p/f(Re) is also exported: this is an instantaneous drag/slip coefficient timescale, not a linear perturbation response time. The longitudinal linearized drag response would involve the derivative of slip*f. No arbitrary gas disturbance duration or atomiser-frequency timescale is assigned to the active Stokes-number results.

For an explicitly hypothetical jet-interaction comparison, the separate ZNMF table supplies its measured local Uc(s) for St_D(s)=tau_p*Uc(s)/0.002. This is a local nozzle-length-based comparison, not nozzle-exit Stokes number, and not a measured PZT Stokes number. Its own s origin is retained without registering it to PZT y. No combined-flow or entrainment claim is made and no ZNMF velocity is substituted into PZT slip. The 190 kHz atomiser frequency is not used as a gas-flow frequency. No frequency-based Stokes number is computed.

## Results and limits

The measured maximum remains 1.686548142 m/s at x=14.54665 mm, y=1.284 mm. The quiescent-air model curves describe conditional deceleration and rise distances for the selected diameters, not fitted or validated experimental trajectories. Diameter, coflow and spatial sampling can all alter comparison; disagreement cannot uniquely identify any one cause. Without measured diameter and simultaneous gas velocity there is no uniquely measured Re_d or PZT Stokes number. The separate ZNMF profile exceeding some droplet velocities does not establish successful entrainment.

## Files and regeneration

Every plotted line has CSV values: measured_fixed_x_profile.csv, ballistic_reference.csv, trajectory_quiescent_d*um.csv, reynolds_vs_height.csv, particle_relaxation_time.csv, diameter_sensitivity.csv and znmf_local_nozzle_stokes.csv. parameters_used.json records properties, assumptions, sources, initial state and calibration provenance. validation.json records signs, units/scale checks, low-Re limit, numerical step convergence and unchanged experimental hashes. Configuration: ../../../config/theory_comparison.json. Reusable source: ../../../analysis/theory_comparison.py.

From workspace root: `.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.theory_comparison`.

Optional simultaneous gas data must use identical native y nodes and documented registration in the configuration. Its two components can produce measured-gas-slip Re sensitivity cases. A spatial trajectory comparison in a nonuniform measured gas field additionally needs a spatially resolved gas interpolant valid along the droplet path; the fixed-x gas column alone is not used to invent one.

## Figure captions

velocity_theory_comparison: Signed fixed-x mean axial droplet velocity versus h=y-y0 relative to the first finite paired measurement, gravity-only ballistic reference and assumed quiescent-air spherical-droplet trajectories for configured diameters. Model branches terminate at their turning points or the upper domain; the ensemble profile is not a trajectory.

reynolds_vs_height: In-plane-slip Reynolds estimates versus h=y-y0 from the measured mean droplet vector under assumed stationary gas and configured diameters. No measured diameter or gas velocity is implied.

particle_relaxation_time: Stokes-limit particle relaxation time versus assumed constant diameter; finite-Re drag requires the separately documented correction.

diameter_sensitivity: Ascending rise within the measured domain for assumed still-air droplets; endpoints are turning points unless marked domain_top in the CSV. The ballistic reference uses the same measured initial axial velocity.

znmf_local_nozzle_stokes: Hypothetical jet-interaction local nozzle-length Stokes comparison using separately measured ZNMF Uc(s), D0=2 mm, and assumed droplet diameters. No PZT gas registration or coflow is implied.

## Sources

- water_density: https://www.nist.gov/system/files/documents/2017/05/09/NISTIR_7383_2013-04-24_2015-07-13Rev.pdf

- air_density: https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=104&refer=ap

- air_viscosity: https://ntrs.nasa.gov/api/citations/20250005707/downloads/main.pdf?attachment=true

- drag_and_response: https://doc.comsol.com/6.4/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.36.html

## Computed sensitivity summary

| Assumed d (micrometres) | Stokes-limit tau (ms) | Rise from y0 (mm) | Maximum trajectory in-plane Re |
|---|---|---|---|
| 5 | 0.0764563 | 0.119127 | 0.552401 |
| 10 | 0.305825 | 0.457435 | 1.1048 |
| 20 | 1.2233 | 1.69778 | 2.2096 |
| 50 | 7.64563 | 8.55511 | 5.52401 |
| 100 | 30.5825 | 24.3443 | 11.048 |
