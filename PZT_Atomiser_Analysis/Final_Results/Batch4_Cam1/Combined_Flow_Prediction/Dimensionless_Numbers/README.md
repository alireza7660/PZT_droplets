# Geometry-independent dimensionless diagnostics

Nominal model diameter 10 um; tau_Stokes=0.000305825036194 s;
documented ZNMF actuation f=120 Hz; Stk_f=0.0366990043432 (no 2*pi).
Acquisition 10 Hz and atomiser 190 kHz are not used as the jet actuation frequency.

CSV files:
- nominal_droplet_response.csv: diameter, Stokes-limit time, nu_air, D0, f, Stk_f.
- pzt_zero_background_reference_field.csv: signed measured mean u/v, coordinates,
  N and conditional particle Re under an explicit zero-background-gas reference.
- pzt_fixed_line_zero_background_reference.csv: original fixed-line profile,
  native indices and the same reference Re. Ambient PZT air was not measured.
- znmf_local_centreline_numbers.csv: all 192 original rows, measured Uc, support
  flags, local Re_D,c, Sr_D,c, Stk_D and Stk_f. Ten unsupported rows remain gaps.

PDF/600-dpi PNG stems:
pzt_Re_p_zero_background_reference; pzt_profile_Re_p_zero_background_reference;
znmf_Re_D_c; znmf_Sr_D_c; znmf_Stk_D; droplet_Stk_f.

Re_D,c=Uc*D0/nu; Sr_D,c=f*D0/Uc; Stk_D=tau_Stokes*Uc/D0. These are local
centreline-based diagnostics, never Re0/Sr0 or measured combined-flow trajectories.
Stokes Stk is distinct from Strouhal Sr. See the single methods_results.tex section,
parameters_and_validation.json and source_manifest.json for values and provenance.

Run independently of ALL placement settings:
.\.venv\Scripts\python -m PZT_Atomiser_Analysis.analysis.prospective_studies --dimensionless-only
