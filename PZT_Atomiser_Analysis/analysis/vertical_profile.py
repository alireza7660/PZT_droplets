"""Native fixed-x axial profile through the finite global mean-v maximum."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .plotting import contour_map,save_figure


STEMS=('mean_v_vertical_profile_location','mean_v_vertical_profile')


def extract_vertical_profile(statistics_path):
    """Return every row of one native column, ordered by increasing y.

    Missing mean velocities remain NaN and are never removed or filled. Ties
    for the finite global maximum select the first native row-major location.
    """
    with np.load(statistics_path) as z:
        x,y,v,n=(z[k].copy() for k in ['x_mm','y_mm','mean_v','N'])
    if not (x.shape==y.shape==v.shape==n.shape) or v.ndim!=2:
        raise ValueError('Expected aligned two-dimensional native fields')
    good=np.isfinite(v)
    if not good.any():raise ValueError('No finite mean axial velocity')
    i,j=np.unravel_index(np.argmax(np.where(good,v,-np.inf)),v.shape)
    if not np.isfinite(y[:,j]).all():raise ValueError('Selected column has undefined y coordinates')
    if not np.allclose(x[:,j],x[i,j],rtol=0,atol=1e-12):
        raise ValueError('Native column is not a fixed-x profile in this coordinate mapping')
    rows=np.argsort(y[:,j],kind='stable')
    profile=dict(y_mm=y[rows,j],x_mm=x[rows,j],mean_v=v[rows,j],N=n[rows,j],
                 native_row=rows.astype(np.int64),native_column=np.full(len(rows),j,dtype=np.int64))
    maximum=dict(native_row=int(i),native_column=int(j),x_mm=float(x[i,j]),
                 y_mm=float(y[i,j]),mean_v_m_s=float(v[i,j]),N=int(n[i,j]),
                 finite_maximum_tie_count=int(np.count_nonzero(good & (v==v[i,j]))),
                 measured_y_range_mm=[float(np.nanmin(y)),float(np.nanmax(y))],
                 profile_y_range_mm=[float(profile['y_mm'][0]),float(profile['y_mm'][-1])])
    assert np.all(np.diff(profile['y_mm'])>=0)
    assert len(profile['mean_v'])==v.shape[0]
    at_peak=np.flatnonzero(profile['native_row']==i)
    assert len(at_peak)==1 and profile['mean_v'][at_peak[0]]==v[i,j]
    np.testing.assert_array_equal(profile['mean_v'],v[profile['native_row'],j])
    return profile,maximum


def profile_methods(maximum):
    return (f'## Fixed-x vertical axial profile\n\nThe finite global maximum of mean_v selects native '
            f'column {maximum["native_column"]} at x={maximum["x_mm"]:.8g} mm. '
            f'The maximum is {maximum["mean_v_m_s"]:.9g} m/s at '
            f'y={maximum["y_mm"]:.8g} mm, native row {maximum["native_row"]}. '
            'Every measured row in that column is retained and sorted together with its y coordinate, '
            'velocity, sample count and native indices. Signed values and gaps are preserved. '
            'No adjacent-column averaging, smoothing, interpolation or extrapolation is applied. '
            'This is a fixed-x axial profile through one global maximum, not a tracked centreline; '
            'the plume may shift laterally with height. The location map retains the existing '
            'mean-v colour limits, levels, linear normalization and physical aspect ratio. '
            'Its thin dashed line spans the measured y domain and its small black marker shows the maximum.\n')


def captions(maximum):
    return {
        STEMS[0]:(rf'Mean axial velocity $\overline{{v}}$ with a thin black dashed line at '
                  rf'$x={maximum["x_mm"]:.5f}\,\mathrm{{mm}}$, spanning the measured domain. '
                  'The black marker identifies the finite global mean-field maximum. The original '
                  'symmetric linear colour scale and all finite signed data are retained.'),
        STEMS[1]:(rf'Fixed-$x$ axial profile of $\overline{{v}}$ at '
                  rf'$x={maximum["x_mm"]:.5f}\,\mathrm{{mm}}$, through the global mean-field maximum. '
                  'One native column is ordered by increasing y over its full measured range; '
                  'signed values and gaps are retained without interpolation, smoothing or column averaging. '
                  'This profile is not a tracked centreline.')}


def update_documents(out,manifest,maximum):
    from .calibrated import write_json
    method=out/'methods/METHODS.md';text=method.read_text(encoding='utf-8')
    text=text.split('## Fixed-x vertical axial profile')[0].rstrip()+'\n\n'+profile_methods(maximum)
    method.write_text(text,encoding='utf-8')
    md=out/'methods/CAPTIONS.md';text=md.read_text(encoding='utf-8')
    text=text.split('## Fixed-x vertical axial profile')[0].rstrip()+'\n\n## Fixed-x vertical axial profile\n\n'
    for name,caption in captions(maximum).items():text+='**'+name+'**: '+caption+'\n\n'
    md.write_text(text,encoding='utf-8')
    tex=out/'methods/CAPTIONS.tex';text=tex.read_text(encoding='utf-8')
    text=text.split('% Fixed-x vertical axial profile')[0].rstrip()+'\n\n% Fixed-x vertical axial profile\n'
    for name,caption in captions(maximum).items():text+='% '+name+'\n'+r'\caption{'+caption+'}\n\n'
    tex.write_text(text,encoding='utf-8')
    write_json(out/'methods/figure_manifest.json',manifest)
    selection={cat:[e['id'] for e in manifest if e['category']==cat] for cat in ['primary','supporting','supplementary','diagnostic']}
    write_json(out/'methods/report_figure_selection.json',selection)
    summary_path=out/'methods/summary.json';summary=json.loads(summary_path.read_text())
    before=summary['figure_count'];summary['figure_count']=len(manifest)
    summary['publication_figure_count']=sum(e['category']!='diagnostic' for e in manifest)
    summary['vertical_profile_selection']=maximum
    write_json(summary_path,summary)
    readme=out/'README.md'
    text=readme.read_text(encoding='utf-8').replace(f'{before} figures,',f'{len(manifest)} figures,')
    readme.write_text(text,encoding='utf-8')


def generate_vertical_profile(out,cfg=None,exporter=None,update_existing_documents=True):
    """Generate only the two new figures and extraction tables from final data."""
    from .final_analysis import Exporter,AXES,LABELS
    from .calibrated import write_csv,write_json
    out=Path(out).resolve();source=out/'data/statistics.npz'
    source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    profile,maximum=extract_vertical_profile(source)
    with np.load(source) as z:x,y,v=(z[k].copy() for k in ['x_mm','y_mm','mean_v'])
    # Reuse the actual existing mean-map scale and exact contour boundaries.
    with np.load(out/'data/plot_data/mean_v.npz') as z:
        color_limits=z['color_limits'].copy();levels=z['levels'].copy()
        np.testing.assert_array_equal(z['values'],v)
        assert str(z['normalization'])=='linear'
    owned=exporter is None
    if owned:
        exporter=Exporter(out)
        exporter.manifest=json.loads((out/'methods/figure_manifest.json').read_text())
    exporter.manifest=[e for e in exporter.manifest if e['id'] not in STEMS]
    lines=[{'x':maximum['x_mm'],'y_min':maximum['measured_y_range_mm'][0],'y_max':maximum['measured_y_range_mm'][1]}]
    points=[{'x':maximum['x_mm'],'y':maximum['y_mm']}]
    stem,data=exporter.paths(STEMS[0])
    contour_map(x,y,v,label=LABELS['mean_v'],limits=color_limits,levels=levels,
                output=stem,**AXES,signed=True,data_output=data,
                vertical_lines=lines,point_markers=points)
    exporter.record(STEMS[0],captions(maximum)[STEMS[0]])
    exporter.manifest[-1]['source_statistics']='data/statistics.npz'
    fig,ax=plt.subplots(figsize=(6.2,4.4),layout='constrained')
    ax.plot(profile['y_mm'],np.ma.masked_invalid(profile['mean_v']),color='#d62728',linewidth=1.3,
            label=rf'$x={maximum["x_mm"]:.5f}\,\mathrm{{mm}}$')
    ax.set(xlabel=r'$y$ ($\mathrm{mm}$)',ylabel=LABELS['mean_v']);ax.legend()
    payload=dict(profile,coordinate_unit='mm',velocity_unit='m/s',selected_x_mm=maximum['x_mm'],
                 global_maximum_y_mm=maximum['y_mm'],global_maximum_mean_v=maximum['mean_v_m_s'])
    stem,_=exporter.paths(STEMS[1]);profile_npz=out/'data'/f'{STEMS[1]}.npz'
    save_figure(fig,stem,payload,data_output=profile_npz)
    write_csv(out/'data'/f'{STEMS[1]}.csv',
              (dict(y_mm=float(profile['y_mm'][k]),x_mm=float(profile['x_mm'][k]),mean_v=float(profile['mean_v'][k]),
                    N=int(profile['N'][k]),native_row=int(profile['native_row'][k]),native_column=int(profile['native_column'][k])) for k in range(len(profile['y_mm']))))
    exporter.record(STEMS[1],captions(maximum)[STEMS[1]])
    exporter.manifest[-1]['numerical_source']=str(profile_npz.relative_to(out))
    exporter.manifest[-1]['source_statistics']='data/statistics.npz'
    for entry in exporter.manifest[-2:]:entry['native_column']=maximum['native_column'];entry['fixed_x_mm']=maximum['x_mm']
    write_json(out/'methods/vertical_profile_selection.json',maximum)
    # Independent geometry and extraction checks against the saved plot sources.
    with np.load(data) as z:
        np.testing.assert_array_equal(z['values'],v)
        np.testing.assert_array_equal(z['color_limits'],color_limits)
        np.testing.assert_array_equal(z['levels'],levels)
        np.testing.assert_allclose(z['vertical_lines_xy'],[[maximum['x_mm'],y.min(),y.max()]],rtol=0,atol=1e-14)
        np.testing.assert_allclose(z['point_markers_xy'],[[x[maximum['native_row'],maximum['native_column']],y[maximum['native_row'],maximum['native_column']]]])
    assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash
    report={'selection':maximum,'profile_point_count':len(profile['mean_v']),
            'all_native_rows_retained':True,'signed_values_and_gaps_preserved':True,
            'fixed_native_column':True,'global_maximum_included':True,
            'location_line_and_marker_verified':True,'existing_colour_scale_preserved':True,
            'source_statistics_sha256_unchanged':source_hash}
    write_json(out/'methods/vertical_profile_validation.json',report)
    if owned and update_existing_documents:update_documents(out,exporter.manifest,maximum)
    return profile,maximum


def notebook_cells():
    import nbformat as nb
    cells=[nb.v4.new_markdown_cell('''## Fixed-x vertical axial profile through the global maximum

The profile uses one native column containing the finite global mean-v
maximum. It is sorted by increasing y, retaining signed values, gaps and local
sample counts. It is not a tracked centreline; no maximum is selected at each
height. Full regeneration produces these figures through the same pipeline.
'''),nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.vertical_profile import extract_vertical_profile
vertical_profile, vertical_maximum=extract_vertical_profile(output/'data/statistics.npz')
print(json.dumps(vertical_maximum,indent=2))
for name in ['mean_v_vertical_profile_location','mean_v_vertical_profile']:
    display(Image(filename=str(output/'figures'/(name+'.png')),width=600))
''')]
    for cell in cells:cell.metadata['fixed_x_vertical_profile']=True
    return cells


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',default='PZT_Atomiser_Analysis/config/analysis.json')
    args=parser.parse_args();config_path=Path(args.config).resolve();cfg=json.loads(config_path.read_text())
    profile,maximum=generate_vertical_profile(config_path.parents[2]/cfg['output'],cfg)
    print(json.dumps(maximum,indent=2))
