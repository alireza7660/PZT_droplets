"""Concise scientific text and exact report-file selection for theory_package."""
import json


def write_report(out,params,manifest):
    cfg=params['configuration'];p=params['properties'];s=params['initial_condition']
    descriptions={
      'velocity_reference_full': (
        r'Fixed-line ensemble-mean axial droplet velocity (red) and individual spherical-droplet reference trajectories in assumed stationary air, with assumed diameters $d=5$, 10, 20, 50 and $100\,\mu\mathrm{m}$. '
        r'The dashed curve is the gravity-only reference. Height $h=y-y_0$ starts at the first measurement used for initialisation, not the atomiser surface. '
        r'The ensemble mean is neither a particle trajectory nor a tracked plume centreline. Reference curves stop at their upward turning point or the measurement boundary; no diameter is inferred from agreement.'),
      'velocity_reference_near_origin': (
        r'Near-origin view, $0\leq h\leq2\,\mathrm{mm}$, of exactly the same data and reference curves, colours and line styles as the full-range velocity figure. '
        r'The small-droplet stationary-air branches are resolved without altering their measured initial velocity, extending stopped droplets, or adding descending branches. Assumed diameters are sensitivity cases.'),
      'reynolds_measured_velocity_assumed_air': (
        r'Conditional in-plane Reynolds estimates $Re_d^{\mathrm{mean}}=\rho_a d\sqrt{\overline{u}^{\,2}+\overline{v}^{\,2}}/\mu_a$ using the measured ensemble-mean droplet vector, assumed stationary air and assumed diameters. '
        r'This is Reynolds number based on the magnitude of the mean velocity vector, not an ensemble average of instantaneous Reynolds numbers. '
        r'Curves span the measured domain because their velocities come from the spatial measurements, independently of the trajectory reference endpoints.'),
      'reynolds_reference_trajectories': (
        r'Reynolds number $Re_d^{\mathrm{ref}}=\rho_a d\lVert\boldsymbol{w}_a-\boldsymbol{w}_d\rVert/\mu_a$ along each individual reference trajectory in assumed stationary air. '
        r'Each assumed-diameter curve ends at its own upward turning point or the measurement boundary. These trajectory values are distinct from the conditional estimates based on measured ensemble means; the curves are not extended across the remaining measurement domain.'),
      'particle_relaxation_time': (
        r'Stokes-limit response time $\tau_{p,\mathrm{Stokes}}=\rho_w d^2/(18\mu_a)$ versus assumed droplet diameter, using water and air properties at $20\,^{\circ}\mathrm{C}$. '
        r'This is a low-Reynolds-number baseline, not an exact finite-Reynolds-number relaxation time. The separately exported $\tau_{\mathrm{drag,sec}}=\tau_{p,\mathrm{Stokes}}/f(Re_d)$ is a drag/slip coefficient timescale.'),
      'diameter_sensitivity': (
        r'Actual upward turning heights of the stationary-air droplet references and the gravity-only reference, measured from the initialisation sample. '
        r'The dotted line marks the upper measurement boundary, which is a plotting limit rather than a stopping height. '
        r'The logarithmic ordinate distinguishes the assumed-diameter cases. Actual turning heights and domain-truncated plot endpoints are exported separately.'),
      'Jet_Interaction_Sensitivity/znmf_local_nozzle_stokes': (
        r'Hypothetical jet-interaction sensitivity $St_D(x)=\tau_{p,\mathrm{Stokes}}U_c(x)/D_0$, with $D_0=2\,\mathrm{mm}$ and assumed droplet diameters. '
        r'$U_c(x)$ is from the separately measured ZNMF jet in its own downstream coordinate. This is a local nozzle-length-based response comparison, not a measured PZT-plume or nozzle-exit Stokes number. '
        r'The downstream trend follows the chosen $U_c(x)/D_0$ timescale. No spatial alignment, simultaneous PZT gas measurement or successful entrainment is established; the atomiser frequency is not used as a gas-flow frequency.')}
    methods=rf'''\section{{Droplet velocities and conditional reference comparisons}}
The PZT atomiser operated at 190 kHz and 9 Vpp. Final calibrated correlation-derived
droplet velocities were used unchanged: $5.35\,\mu\mathrm{{m/pixel}}$,
$\Delta t=25\,\mu\mathrm{{s}}$, and $0.214\,\mathrm{{m\,s^{{-1}}}}$ per pixel
displacement, applied once in the experimental analysis. Coordinates $x,u$ are
horizontal/rightward and $y,v$ upward. The measured profile is the ensemble mean
$\overline{{v}}(x_*,y)$ along the existing fixed line $x_*={s['fixed_x_mm']:.8g}\,\mathrm{{mm}}$;
it is neither an individual droplet trajectory nor a tracked plume centreline.
We use $h=y-y_0$, with image-reference $y_0={s['initial_y_mm']:.8g}\,\mathrm{{mm}}$
at the first finite paired measurement used for initialisation. It is not an
established atomiser surface. The initial measured axial velocity is
$v_0={s['initial_v_m_s']:.9g}\,\mathrm{{m\,s^{{-1}}}}$; the measured horizontal
mean is also retained in the two-component reference calculation.

Individual spherical water-droplet references are calculated in assumed stationary
air, with constant assumed diameters 5, 10, 20, 50 and $100\,\mu\mathrm{{m}}$,
neglecting evaporation. The air temperature is $20\,^{{\circ}}\mathrm{{C}}$;
pressure 101.325 kPa and equal water/air temperature are assumptions. Humidity is
unavailable; dry-air properties are an approximation. Adopted values are
$\rho_a={p['rho_air_kg_m3']:.7g}\,\mathrm{{kg\,m^{{-3}}}}$,
$\rho_w={p['rho_water_kg_m3']:.8g}\,\mathrm{{kg\,m^{{-3}}}}$ and
$\mu_a={p['mu_air_Pa_s']*1e5:.7g}\times10^{{-5}}\,\mathrm{{Pa\,s}}$.
Water density follows NIST GLP10, air density uses a NIST dry-air reference,
and viscosity uses the NASA-reported Sutherland constants
\cite{{pzt-water,pzt-air,pzt-viscosity}}.
\begin{{align}}
\dot{{\boldsymbol q}}&=\boldsymbol w_d,\\
\dot{{\boldsymbol w}}_d&=-g(1-\rho_a/\rho_w)\boldsymbol e_y
+\frac{{f(Re_d)}}{{\tau_{{p,\mathrm{{Stokes}}}}}}
(\boldsymbol w_a-\boldsymbol w_d),\qquad \boldsymbol w_a=\boldsymbol 0,\\
Re_d&=\frac{{\rho_a d\lVert\boldsymbol w_a-\boldsymbol w_d\rVert}}{{\mu_a}},
\qquad\tau_{{p,\mathrm{{Stokes}}}}=\frac{{\rho_w d^2}}{{18\mu_a}}.
\end{{align}}
Here $g=9.80665\,\mathrm{{m\,s^{{-2}}}}$ and the spherical drag multiplier is
$f=1$ for $Re_d\leq0.01$ and $f=1+0.15Re_d^{{0.687}}$ otherwise, with
$Re_d<800$ enforced \cite{{pzt-drag}}. Only in-plane slip is available. Rigid-interface spherical
drag, fixed diameter and one-way coupling are model assumptions; deformation,
internal circulation, evaporation and particle interactions are omitted.
The gravity-only reference omits drag and buoyancy. Each displayed upward
trajectory terminates at its turning point or the measurement boundary.

For assumed diameters 5, 10, 20, 50 and $100\,\mu\mathrm{{m}}$, the stationary-air
turning heights are '''
    stops=[e['actual_turning_height_mm'] for e in s['endpoints'][1:]]
    methods+=', '.join(f'{h:.5g}' for h in stops)+rf'''$\,\mathrm{{mm}}$, respectively.
The gravity-only turning height is ${s['endpoints'][0]['actual_turning_height_mm']:.6g}\,\mathrm{{mm}}$,
whereas the plotted measurement extent is only ${s['measurement_height_extent_mm']:.6g}\,\mathrm{{mm}}$.
These are different quantities. The full and $0$--$2\,\mathrm{{mm}}$ views
show the same unsmoothed measurements and unchanged reference curves.
Disagreement with stationary-air trajectories does not itself identify a coding
error or determine droplet size. Possible contributors include surrounding gas
motion, size and initial-velocity distributions, and lateral plume displacement;
the available evidence does not establish which dominates. Neither diameter nor
an unknown gas field is fitted, and the reference model is not experimentally validated.

Two Reynolds-number quantities are reported separately. The conditional estimate
$Re_d^{{\mathrm{{mean}}}}=\rho_a d\lVert\overline{{\boldsymbol w}}_d\rVert/\mu_a$
uses the measured mean vector and assumed stationary air; it is not
$\overline{{Re_d}}$ from instantaneous samples and spans the measurement domain.
The trajectory reference $Re_d^{{\mathrm{{ref}}}}$ uses the model's instantaneous
velocity and ends at that trajectory's endpoint. They are not equivalent results.
The response-time plot reports the Stokes-limit $\tau_{{p,\mathrm{{Stokes}}}}$.
At finite Reynolds number the exported $\tau_{{\mathrm{{drag,sec}}}}=
\tau_{{p,\mathrm{{Stokes}}}}/f(Re_d)$ is an instantaneous drag/slip coefficient
timescale, not an exact or linearised relaxation time.

The separate hypothetical jet-interaction comparison uses
$St_D(x)=\tau_{{p,\mathrm{{Stokes}}}}U_c(x)/D_0$, $D_0=2\,\mathrm{{mm}}$,
where $U_c(x)$ is from the separately measured ZNMF jet and its own downstream
coordinate. The assumed diameters and the local timescale $D_0/U_c(x)$ define a
local nozzle-length-based response comparison; the downstream trend reflects
that chosen timescale. This is neither a measured PZT-plume Stokes number nor
a nozzle-exit Stokes number. No spatial registration to the PZT field is assumed,
the separate jet velocity is not inserted into the PZT droplet slip, and the
calculation does not demonstrate entrainment. The 190 kHz atomiser frequency is
not a gas-flow frequency in this analysis.
'''
    (out/'methods_results.tex').write_text(methods,encoding='utf-8')
    bibliography=r'''@techreport{pzt-water,
  author = {{National Institute of Standards and Technology}},
  title = {Selected Procedures for Volumetric Calibrations: GLP 10 Water Density Calculations},
  number = {NISTIR 7383}, year = {2013},
  url = {https://www.nist.gov/system/files/documents/2017/05/09/NISTIR_7383_2013-04-24_2015-07-13Rev.pdf}
}
@misc{pzt-air,
  author = {{National Institute of Standards and Technology}},
  title = {Composition of Air, Dry (Near Sea Level)},
  url = {https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=104&refer=ap}
}
@misc{pzt-viscosity,
  author = {{NASA}}, title = {NASA Technical Report 20250005707},
  note = {Sutherland viscosity constants: mu0=1.716e-5 Pa s, T0=273.15 K, S=110.4 K},
  url = {https://ntrs.nasa.gov/api/citations/20250005707/downloads/main.pdf?attachment=true}
}
@misc{pzt-drag,
  author = {{COMSOL}}, title = {Particle Motion in a Fluid},
  note = {COMSOL 6.4 documentation, Stokes and Schiller--Naumann drag laws},
  url = {https://doc.comsol.com/6.4/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.36.html}
}
'''
    (out/'references.bib').write_text(bibliography,encoding='utf-8')
    captions='';environments='% Requires graphicx. Compile from Theory_Comparison or set graphicspath accordingly.\n'
    for item in manifest:
        stem=item['stem'];caption=descriptions[stem]
        captions+='% '+stem+'\n'+r'\caption{'+caption+'}\n\n'
        environments+=r'\begin{figure}[htbp]\centering'+'\n'+r'\includegraphics[width=0.92\linewidth]{'+stem+'.pdf}\n'+r'\caption{'+caption+'}\n'+r'\label{fig:'+stem.replace('/','-').replace('_','-')+'}\n'+r'\end{figure}'+'\n\n'
    (out/'captions.tex').write_text(captions,encoding='utf-8')
    (out/'report_figures.tex').write_text(environments,encoding='utf-8')
    selection=dict(figures=[dict(stem=f['stem'],pdf=f['stem']+'.pdf',png=f['stem']+'.png',category=f['category']) for f in manifest],
                   text=['methods_results.tex','captions.tex','report_figures.tex','references.bib'],
                   excluded=['Archive_Before_Report_Revision_1522f54/','validation/'],
                   notes='Only listed active files are selected for the scientific report; archives and computational checks are excluded.')
    (out/'report_manifest.json').write_text(json.dumps(selection,indent=2))
    readme=f'''# Fixed-line droplet mean velocities and stationary-air references

This package distinguishes a measured ensemble-mean profile from individual
monodisperse reference trajectories. Assumed diameters are {cfg['diameters_um']}
micrometres. They are not measured or fitted sizes. Air at 20 C is recorded;
stationary air, pressure 101325 Pa and equal water temperature are assumptions.
Humidity is unavailable; dry-air properties are approximate and evaporation is omitted.

## Measurement and references

The native fixed line is x={s['fixed_x_mm']:.8g} mm. It is not a tracked centreline.
Height h=y-y0 starts at y0={s['initial_y_mm']:.8g} mm (image reference), with measured
initial v={s['initial_v_m_s']:.9g} m/s and u={s['initial_u_m_s']:.9g} m/s. The location
is selected by finite paired means and N>={cfg['minimum_initial_sample_count']}; this
is sample support, not independent physical validation or an atomiser-surface location.
Every profile row, signed value, gap, native index and original y coordinate is retained.
Existing final calibration and sample selection are unchanged; SI velocities are not recalibrated.

Gravity, buoyancy and spherical finite-Re drag define the stationary-air references.
The same physical equations, properties and measured initial state used before this
revision are retained. Each upward branch stops at v=0 or the measured upper boundary.
`reference_endpoints.csv` separates actual turning heights from plotted endpoints.
Gravity-only actual rise is {s['endpoints'][0]['actual_turning_height_mm']:.8g} mm;
its plotted curve stops at the {s['measurement_height_extent_mm']:.8g} mm measurement
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

'''
    for f in manifest:readme+=f"- `{f['stem']}.pdf` (also 600-dpi PNG and SVG).\n"
    readme+='''
Every curve has CSV provenance in `figure_curve_manifest.json`; figures are drawn
directly from those CSV values. Per-diameter trajectory CSVs include original y,
relative height h, time, vector components, in-plane Re and defined timescales.
`diameter_sensitivity.csv` reports actual versus displayed endpoints; the full
endpoint table also includes gravity-only motion. `parameters_used.json` records
all assumptions, properties, source URLs, calibration, initial state and ZNMF provenance.

## Regeneration

From the PZT_droplets workspace root, run this single command using the project environment:

```powershell
.\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.update_theory_notebook
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
'''
    (out/'README.md').write_text(readme,encoding='utf-8')
    jet_readme=r'''# Hypothetical jet-interaction sensitivity

St_D(x)=tau_p,Stokes*Uc(x)/D0, D0=0.002 m.
tau_p,Stokes=rho_water*d^2/(18*mu_air) is the Stokes-limit response-time baseline;
the assumed diameters are 5, 10, 20, 50 and 100 micrometres. It is not an exact
finite-Re relaxation time. Uc(x) is from the separately measured ZNMF jet; its
own downstream coordinate is retained without physical registration to PZT y.
The local gas timescale is D0/Uc(x). Its variation determines the downstream
trend in this local nozzle-length-based comparison. This is neither a measured
PZT-plume Stokes number nor a nozzle-exit Stokes number. It does not demonstrate
successful entrainment or supply gas velocity measured in the droplet experiment.
No 190 kHz gas-flow frequency is assumed. The logarithmic ordinate resolves all
assumed-diameter cases. See parameters_used.json and the plotted-value CSV here.
Regenerate using the single command in the parent README. No separate physical
alignment or hypothetical combined-flow trajectory is inferred.
'''
    (out/'Jet_Interaction_Sensitivity/README.md').write_text(jet_readme,encoding='utf-8')
