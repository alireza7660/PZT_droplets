"""Mask-aware, unsmoothed PIV maps; shared limits are set before plotting."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.linewidth': .8, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'savefig.dpi': 600})


def shared_limits(fields, signed=False, correlation=False):
    """Use full finite range, never percentile clipping, across all conditions."""
    if correlation:
        return (-1., 1.)
    values = [np.asarray(z)[np.isfinite(z)] for z in fields]
    values = [z for z in values if z.size]
    if not values:
        raise ValueError('No finite data for color limits')
    lo, hi = min(z.min() for z in values), max(z.max() for z in values)
    if signed:
        extent = max(abs(lo), abs(hi), np.finfo(float).eps)
        return (-float(extent), float(extent))
    return (0., float(max(hi, np.finfo(float).eps)))


def save_figure(fig, stem, data):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(stem.with_suffix('.' + ext), bbox_inches='tight', dpi=600)
    np.savez_compressed(stem.with_suffix('.npz'), **data)
    plt.close(fig)


def contour_map(x, y, values, *, label, limits, output, xlabel, ylabel,
                signed=False, levels=41):
    """Piecewise contour interpolation only within fully valid native grid cells.

    No filtering, filling, extrapolation or triangulation across missing nodes.
    corner_mask=False removes any quad touching an invalid node. Contour level
    interpolation is intrinsic to contourf, not a spatial smoothing operation.
    """
    z = np.ma.masked_invalid(np.ma.asarray(values, dtype=float))
    fig, ax = plt.subplots(figsize=(6.2, 5.0), layout='constrained')
    m = ax.contourf(x, y, z, levels=np.linspace(*limits, levels),
                    cmap='RdBu_r' if signed else 'viridis',
                    corner_mask=False, antialiased=False)
    fig.colorbar(m, ax=ax, label=label)
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect='equal')
    save_figure(fig, output, {'x': x, 'y': y, 'values': z.filled(np.nan),
                             'color_limits': limits, 'levels': np.linspace(*limits, levels)})


def quiver_map(x, y, transverse, axial, *, output, xlabel, ylabel,
               velocity_unit, reference_speed, stride=3):
    """Separate vectors on the measured grid, with a common physical arrow scale."""
    valid = np.isfinite(transverse) & np.isfinite(axial)
    u = np.ma.array(transverse, mask=~valid)
    v = np.ma.array(axial, mask=~valid)
    sl = (slice(None, None, stride), slice(None, None, stride))
    fig, ax = plt.subplots(figsize=(6.2, 5.0), layout='constrained')
    width = max(np.ptp(x), np.ptp(y))
    q = ax.quiver(x[sl], y[sl], u[sl], v[sl], angles='xy',
                  scale_units='xy', scale=reference_speed/(.07*width), width=.003)
    key = ax.quiverkey(q, .72, .93, reference_speed,
                      f'{reference_speed:.3g} {velocity_unit}', labelpos='S',
                      coordinates='axes')
    key.text.set_bbox(dict(facecolor='white',edgecolor='none',alpha=.9,pad=2))
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect='equal',
           xlim=(np.nanmin(x), np.nanmax(x)), ylim=(np.nanmin(y), np.nanmax(y)))
    save_figure(fig, output, {'x': x, 'y': y, 'transverse': u.filled(np.nan),
                             'axial': v.filled(np.nan), 'stride': stride,
                             'reference_speed': reference_speed})
