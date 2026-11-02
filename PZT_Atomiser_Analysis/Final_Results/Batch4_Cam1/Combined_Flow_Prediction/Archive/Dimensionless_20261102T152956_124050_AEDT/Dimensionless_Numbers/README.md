# Dimensionless-number audit

| Quantity | Value/range |
|---|---|
| Re_D,c | 11.818756–356.18752 |
| Sr_D,c | 0.089536257–2.698397 |
| tau_p (ms) | 0.3058250362 |
| Stk_f | 0.03669900434 |
| Stk_D,c | 0.013600298–0.4098787 |
| Re_p | 2.7042404e-05–0.5125812 |
| Stk_D,path | 0.11698913–0.39034595 |
| Re0, St0 | Unresolved |
| PZT-only measured Re_p | Not identifiable |

Ranges span the full valid measured jet centreline or all four complete supported paths, as applicable. No placement or trajectory was recalculated or changed. h=0 is the first supported measured plume-axis plane, placed at Xj=9.995298977 mm where interpolated Uc=0.960000000 m/s. The sampled centreline field is not rescaled to match that curve exactly.

At 20 C: rho_air=1.20479 kg/m3; mu_air=1.81332212036e-05 Pa s; nu_air=1.50509393368e-05 m2/s; rho_water=998.20674556 kg/m3. Atmospheric pressure 101325 Pa and dry-air properties are assumptions; humidity was not supplied. d=10 micrometres is a nominal model estimate: 1.89 times an assumed 5 micrometre emerging jet gives 9.45 micrometres, rounded to 10; neither jet/aperture equality nor measured droplet sizing is established.

The full valid centreline maximum is 2.68047838041 m/s at x/D0=25.2096, Xj=50.4192 mm. Coordinates x/D0 are multiplied once by D0 in mm. Forcing is 120 Hz; acquisition 10 Hz is never used for Sr or Stk_f. Nonpositive velocities and source gaps remain undefined in local centreline diagnostics. No signed velocity is replaced by its magnitude. Norms occur only in vector-magnitude definitions of Re_p and Stk_D,path.

Stk_f uses forcing period 1/f; Stk_D,c uses centreline transit time D0/Uc; Stk_D,path uses local gas transit time D0/|u_g|. They answer different response questions. All use Stokes-limit tau_p. The solver retains finite-Re factor F; tau_eff=tau_p/F is a separately defined coefficient time, not an exact constant finite-Re relaxation time. Stk<1 does not prove entrainment, capture or successful delivery. The prescribed steady measured mean gas is one-way coupled; predicted paths are not measured combined-flow trajectories. PZT-only droplet velocities do not establish surrounding gas velocity or a measured slip Reynolds number.

Nozzle-reference finding: No documented phase-resolved nozzle-exit velocity over a forcing period or defensible expulsion-averaged reference speed in available experimental metadata. See nozzle_reference_audit.json for exact evidence and paper definition. Downstream maximum, entrance speed and signed full-cycle mean are not substituted for U0.

Reference-station values are in data/reference_stations.csv. Comparison endpoint uses component-wise interpolation on the saved path, followed by recalculation of norms; the refined minimum uses the existing supported event state. Path tables include every saved integration point through supported termination.

Reproduce without rerunning trajectories:
`.\.venv\Scripts\python.exe -m PZT_Atomiser_Analysis.analysis.dimensionless_audit --notebook`
