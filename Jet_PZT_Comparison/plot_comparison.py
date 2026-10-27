"""Overlay final numerical profiles; Python 3, numpy and matplotlib required."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def read(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--znmf', required=True, type=Path)
    p.add_argument('--pzt', required=True, type=Path)
    p.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    args = p.parse_args()
    paths = [args.znmf.resolve(), args.pzt.resolve()]
    hashes = [hashlib.sha256(s.read_bytes()).hexdigest() for s in paths]
    z, t = map(read, paths)
    # Units are determined from explicit column names, never coordinate magnitude.
    keys = [k for k in ('s_mm', 'x_mm', 'x_over_D0') if k in z[0]]
    if len(keys) != 1:
        raise ValueError('ZNMF must have exactly one coordinate column: s_mm, x_mm or x_over_D0')
    key = keys[0]
    factor = 2.0 if key == 'x_over_D0' else 1.0
    zx = np.array([float(r[key]) for r in z]) * factor
    zu = np.array([float(r['Uc_m_per_s']) for r in z])
    measured = np.array([r['status'] == 'measured' for r in z])
    peak = np.nanargmax(np.where(measured, zu, np.nan))
    np.testing.assert_allclose([zx[peak]/2, zu[peak]], [25.2096, 2.680478380412314], rtol=0, atol=1e-8)
    tx = np.array([float(r['y_mm']) for r in t])
    tu = np.array([float(r['mean_v']) for r in t])
    np.testing.assert_allclose([float(r['x_mm']) for r in t], 14.54665, rtol=0, atol=1e-10)
    zm = measured & np.isfinite(zu) & (zx >= 0) & (zx <= 27)
    tm = (tx >= 0) & (tx <= 27)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for name, rows, x, u, mask in [('znmf', z, zx, zu, zm), ('pzt', t, tx, tu, tm)]:
        indices = np.flatnonzero(mask)
        assert np.all(np.diff(x[mask]) > 0)
        dest = out / (name + '_plotted_data.csv')
        fields = ['source_row_index', 's_mm', 'mean_axial_velocity_m_s'] + list(rows[0])
        with dest.open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for i in indices:
                w.writerow(dict(source_row_index=int(i), s_mm=repr(float(x[i])),
                                mean_axial_velocity_m_s=repr(float(u[i])), **rows[i]))
        saved = read(dest)
        np.testing.assert_array_equal([float(r['s_mm']) for r in saved], x[mask])
        np.testing.assert_array_equal([float(r['mean_axial_velocity_m_s']) for r in saved], u[mask])
        records.append(dict(dataset=name, count=len(indices), s_range_mm=[float(x[mask].min()), float(x[mask].max())],
                            maximum_m_s=float(np.nanmax(u[mask])), csv_roundtrip_exact=True))
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 11,
                         'mathtext.fontset': 'dejavuserif', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(6.4, 4.4), layout='constrained')
    ax.plot(zx[zm], zu[zm], color='red', lw=1.5, label='ZNMF jet')
    ax.plot(tx[tm], np.ma.masked_invalid(tu[tm]), color='blue', ls='--', lw=1.5, label='PZT droplets')
    top = np.ceil(max(np.nanmax(zu[zm]), np.nanmax(tu[tm])) * 1.10 / 0.2) * 0.2
    ax.set(xlim=(0, 27), ylim=(0, top), xlabel=r'Axial coordinate, $s$ (mm)',
           ylabel=r'Mean axial velocity (m s$^{-1}$)')
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=2.5)
    ax.legend(frameon=False, loc='upper right')
    assert tuple(ax.get_xlim()) == (0, 27)
    for line, x, u in zip(ax.lines, [zx[zm], tx[tm]], [zu[zm], tu[tm]]):
        np.testing.assert_array_equal(line.get_xdata(), x)
        np.testing.assert_array_equal(np.asarray(line.get_ydata()), u)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(out / ('mean_axial_velocity_comparison.' + ext), dpi=600)
    plt.close(fig)
    assert hashes == [hashlib.sha256(s.read_bytes()).hexdigest() for s in paths]
    report = dict(inputs=[dict(path=str(s), sha256=h) for s,h in zip(paths, hashes)],
                  znmf_coordinate_column=key, znmf_coordinate_multiplier_to_mm=factor,
                  znmf_full_domain_maximum=dict(x_over_D0=float(zx[peak]/2), velocity_m_s=float(zu[peak])),
                  axis_limits_mm=[0,27], vertical_limits_m_s=[0,float(top)], datasets=records,
                  plotted_arrays_match_sources=True, original_sources_unchanged=True,
                  velocity_recalibration_applied=False)
    (out/'verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
