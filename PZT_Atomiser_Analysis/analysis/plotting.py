"""Mask-aware, unsmoothed PIV maps; shared limits are set before plotting."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.linewidth': .8, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'savefig.dpi': 600,
                     'xtick.direction':'in','ytick.direction':'in',
                     'xtick.top':True,'ytick.right':True,
                     'axes.spines.top':True,'axes.spines.right':True})


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


def save_figure(fig, stem, data, *, data_output=None):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(stem.with_suffix('.' + ext), bbox_inches='tight', dpi=600)
    target=Path(data_output) if data_output is not None else stem.with_suffix('.npz')
    target.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(target, **data)
    plt.close(fig)


def contour_boundaries(limits,levels=41,signed=False):
    """Uniform linear levels; symmetric signed ranges include zero exactly."""
    if np.isscalar(levels):
        count=max(3,int(levels))
        if signed and count%2==0:count+=1
        boundaries=np.linspace(*limits,count)
        if signed:
            if not np.isclose(limits[0],-limits[1]):raise ValueError('Signed linear colour limits must be symmetric about zero')
            boundaries[count//2]=0.
        return boundaries
    boundaries=np.asarray(levels,float)
    np.testing.assert_allclose(np.diff(boundaries),np.diff(boundaries)[0],rtol=1e-12,atol=1e-14)
    return boundaries


def contour_map(x, y, values, *, label, limits, output, xlabel, ylabel,
                signed=False, levels=41, data_output=None,
                markers=None, vertical_lines=None, point_markers=None):
    """Piecewise contour interpolation only within fully valid native grid cells.

    No filtering, filling, extrapolation or triangulation across missing nodes.
    corner_mask=False removes any quad touching an invalid node. Contour level
    interpolation is intrinsic to contourf, not a spatial smoothing operation.
    """
    raw=np.ma.masked_invalid(np.ma.asarray(values,dtype=float))
    display_mask=np.ma.getmaskarray(raw).copy()
    z=np.ma.array(raw.data,mask=display_mask)
    boundaries=contour_boundaries(limits,levels,signed)
    fig, ax = plt.subplots(figsize=(6.2, 5.0), layout='constrained')
    norm=Normalize(vmin=limits[0],vmax=limits[1])
    m = ax.contourf(x, y, z, levels=boundaries,norm=norm,
                    cmap='RdBu_r' if signed else 'viridis',
                    corner_mask=False, antialiased=False)
    cb=fig.colorbar(m, ax=ax, label=label)
    if signed and limits[0]<0<limits[1]:
        cb.set_ticks(np.linspace(limits[0],limits[1],7))
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect='equal')
    for marker in markers or []:
        ax.plot(marker['x'],marker['y'],'o',mfc='white',mec='black',ms=4)
        ax.annotate(marker['id'],(marker['x'],marker['y']),xytext=(5,4),textcoords='offset points',fontsize=8)
    for line in vertical_lines or []:
        ax.plot([line['x'],line['x']],[line['y_min'],line['y_max']],
                color='black',linestyle='--',linewidth=.8)
    for marker in point_markers or []:
        ax.plot(marker['x'],marker['y'],'o',color='black',markersize=3.5)
    save_figure(fig, output, {'x': x, 'y': y, 'values': raw.filled(np.nan),
                             'display_mask':display_mask,'finite_values_only':True,
                             'color_limits':limits,'levels':boundaries,
                             'zero_centred_diverging':signed and limits[0]<0<limits[1],
                             'normalization':'linear','symmetric_signed_limits':signed,
                             'uniform_contour_levels':True,
                             'vertical_lines_xy':np.asarray([[a['x'],a['y_min'],a['y_max']] for a in vertical_lines or []],float).reshape(-1,3),
                             'point_markers_xy':np.asarray([[a['x'],a['y']] for a in point_markers or []],float).reshape(-1,2)},data_output=data_output)


def quiver_map(x, y, transverse, axial, *, output, xlabel, ylabel,
               velocity_unit, reference_speed, stride=3, data_output=None,
               component_names=('transverse','axial')):
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
    save_figure(fig, output, {'x': x, 'y': y, component_names[0]: u.filled(np.nan),
                             component_names[1]: v.filled(np.nan), 'stride': stride,
                             'reference_speed': reference_speed},data_output=data_output)
