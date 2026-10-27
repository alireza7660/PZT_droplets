"""One-time terminology migration; no component or coordinate definitions change."""
from pathlib import Path
import json
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'provisional_calibrated'
    # Preserve the exact numerical baseline and old figure stems for cleanup.
    baseline=out/'diagnostics/terminology_baseline'
    baseline.mkdir(exist_ok=True)
    for p in (out/'results').glob('*_statistics.npz'):
        with np.load(p) as z:
            np.savez_compressed(baseline/p.name,**{k:z[k] for k in z.files})
    prior=json.loads((out/'results/figure_manifest.json').read_text())
    (baseline/'prior_figure_manifest.json').write_text(json.dumps(prior,indent=2))
    p=root/'analysis/calibrated.py'
    s=p.read_text(encoding='utf-8').replace('speed','velocity_magnitude')
    s=s.replace('Mean instantaneous velocity_magnitude (m/s)','Mean instantaneous velocity magnitude (m/s)')
    s=s.replace('Instantaneous velocity_magnitude (m/s)','Instantaneous velocity magnitude (m/s)')
    s=s.replace('In-plane velocity_magnitude (m/s)','In-plane velocity magnitude (m/s)')
    s=s.replace('mean instantaneous velocity_magnitude are distinct','mean instantaneous velocity magnitude are distinct')
    s=s.replace("Transverse position $y$ (mm)","Transverse coordinate, $y$ (mm)")
    s=s.replace("Axial position $x$ (mm)","Axial coordinate, $x$ (mm)")
    start=s.index("METHODS = '''");end=s.index("\ndef run(",start)
    s=s[:start]+s[start:end].replace('velocity_magnitude','velocity magnitude')+s[end:]
    p.write_text(s,encoding='utf-8')
    p=root/'analysis/build_calibrated_notebook.py'
    s=p.read_text(encoding='utf-8').replace('speed','velocity_magnitude')
    s=s.replace('RMS, velocity_magnitude and','RMS, velocity magnitude and')
    p.write_text(s,encoding='utf-8')
    p=root/'README.md'
    s=p.read_text(encoding='utf-8-sig').replace('speed','velocity magnitude')
    s+='\n## Coordinate notation\n\nThe upward spray axis is axial x; horizontal y is transverse, perpendicular to that axis. This assignment follows the experimenter-specified convention, not a universal x/y rule. U is the axial velocity component and V the transverse component. Plots label coordinates separately from velocity components.\n\nVelocity magnitude is the scalar norm sqrt(U squared + V squared). Its ensemble mean differs from the magnitude of the mean velocity vector. Both remain explicitly distinguished.\n'
    p.write_text(s,encoding='utf-8')
    print('Updated active labels, documentation and output names; numerical definitions unchanged.')


if __name__=='__main__':main()
