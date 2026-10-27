"""Redraw signed contours from unchanged final numerical sources, in place."""
import base64
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np
import nbformat as nb
from .final_analysis import LABELS,AXES,UNIT,RENAME,FULL_FIELD_NOTE,limits,documentation,validate
from .plotting import contour_map,contour_boundaries

OLD_NOTE='Negative axial velocities are masked only in this contour display and the colour scale starts at zero. All signed values are retained in the numerical data and statistical calculations; no negative value is replaced by zero or an absolute value.'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    cfg=json.loads((root/'config/analysis.json').read_text());out=root.parent/cfg['output']
    baseline={str(p.relative_to(out)):digest(p) for p in (out/'data').iterdir() if p.is_file()}
    critical=['confirmed_spatial_calibration_um_per_pixel','confirmed_pulse_separation_seconds',
              'coordinate_scale','velocity_scale','components','component_mapping','coordinate_mapping',
              'transform','origin','valid_flag_values','stations','profile_rows','representative_quantiles']
    previous=json.loads((out/'methods/configuration_used.json').read_text())
    assert {k:cfg[k] for k in critical}=={k:previous[k] for k in critical}
    with np.load(out/'data/statistics.npz') as z:
        fields={k:z[k].copy() for k in RENAME};r={'xy':np.stack([z['x_mm'],z['y_mm']]),'N':z['N'].copy(),'frames':{}}
    old_mean=fields['mean_v'].copy()
    summary=json.loads((out/'methods/summary.json').read_text())
    for frame in summary['selected_frame_keys']:
        with np.load(out/'data'/f'frame_{frame}.npz') as z:r['frames'][frame]=np.stack([z['u'],z['v']])
    colour_path=out/'methods/colour_limits.json';colour=json.loads(colour_path.read_text())
    mainlimits=colour['main_ensemble']
    for k in ['mean_v','mean_u','covariance','correlation']:mainlimits[k]=limits([fields[k]],signed=True)
    instant=colour['instantaneous_selected_frames']
    for i,k in enumerate(['u','v']):instant[k]=limits([q[i] for q in r['frames'].values()],signed=True)
    manifest=json.loads((out/'methods/figure_manifest.json').read_text())
    replacements={};changed=[];mask_counts={}
    for entry in manifest:
        name=entry['id'];entry.pop('negative_axial_display_mask',None);entry['display_policy']='all_finite_values'
        entry['caption']=entry['caption'].replace(OLD_NOTE,'').replace(FULL_FIELD_NOTE,'').strip()
        axial=name in ['mean_v','sampling_stations'] or (name.startswith('frame_') and name.endswith('_v')) or (name.startswith('flag_') and name.endswith('_mean_v'))
        signed=axial or name in ['mean_u','covariance','correlation'] or (name.startswith('frame_') and name.endswith('_u')) or (name.startswith('flag_') and name.endswith('_mean_u'))
        if axial:entry['caption']+=' '+FULL_FIELD_NOTE
        if not signed:continue
        source=out/entry['numerical_source'];png=out/entry['files']['png']
        before_png=base64.b64encode(png.read_bytes()).decode('ascii')
        with np.load(source) as z:x,y,values=z['x'].copy(),z['y'].copy(),z['values'].copy()
        markers=None
        if name in fields:lim=mainlimits[name];label=LABELS[name]
        elif name=='sampling_stations':
            lim=mainlimits['mean_v'];label=LABELS['mean_v']
            markers=[{'x':float(x[s['row'],s['column']]),'y':float(y[s['row'],s['column']]),'id':s['id']} for s in cfg['stations']]
        elif name.startswith('frame_'):
            component=name.rsplit('_',1)[1];lim=instant[component];label=rf'${component}$ ({UNIT})'
        else:
            component='mean_v' if name.endswith('_mean_v') else 'mean_u'
            lim=limits([values],signed=True);label=LABELS[component]
        with np.load(source) as z:
            current=(bool(z.get('finite_values_only',False)) and 0. in z['levels']
                     and np.array_equal(z['display_mask'],~np.isfinite(values))
                     and np.array_equal(z['color_limits'],lim)
                     and np.array_equal(z['levels'],contour_boundaries(lim,cfg['contour_levels'],True)))
        if not current:
            contour_map(x,y,values,label=label,limits=lim,output=png.with_suffix(''),**AXES,
                        signed=True,levels=cfg['contour_levels'],data_output=source,markers=markers)
        with np.load(source) as after:
            np.testing.assert_array_equal(after['values'],values)
            np.testing.assert_array_equal(after['x'],x);np.testing.assert_array_equal(after['y'],y)
            np.testing.assert_array_equal(after['display_mask'],~np.isfinite(values))
            assert 0. in after['levels']
            finite=values[np.isfinite(values)]
            assert lim[0]<=finite.min() and lim[1]>=finite.max()
            mask_counts[name]=int(after['display_mask'].sum())
        replacements[before_png]=base64.b64encode(png.read_bytes()).decode('ascii')
        changed.append(name)
    write=lambda path,obj:Path(path).write_text(json.dumps(obj,indent=2),encoding='utf-8')
    write(colour_path,colour);write(out/'methods/configuration_used.json',cfg)
    ex=SimpleNamespace(manifest=manifest)
    documentation(r,fields,cfg,out,ex,mainlimits)
    validate(r,fields,cfg,out,ex)
    # Preserve executed numerical cells; update explanatory cells and actual previews.
    target=root/'notebooks/PZT_final_analysis.ipynb';b=nb.read(target,as_version=4)
    b.cells[0].source='''# PZT atomiser final analysis

Calibration remains 7.25 micrometres/pixel, pulse separation 25 microseconds,
and conversion 0.290 m/s per pixel displacement. Cartesian x/u is horizontal/
transverse, y/v vertical/axial. Existing signs, array orientation, origin and
finite-pair population are unchanged.

All finite signed axial velocities are displayed, including negative values,
on mean, instantaneous and station-location maps. Blue denotes negative,
white zero and red positive. Limits are actual asymmetric data ranges; zero
is an exact contour boundary. Comparable instantaneous maps share limits.
Only genuinely missing or invalid data are masked.
'''
    refreshed=0
    for cell in b.cells:
        if cell.cell_type!='code':continue
        assert cell.execution_count is not None
        assert not any(o.output_type=='error' for o in cell.outputs)
        for obj in cell.outputs:
            if obj.output_type not in ['display_data','execute_result']:continue
            image_data=obj.get('data',{}).get('image/png')
            if image_data in replacements:obj.data['image/png']=replacements[image_data];refreshed+=1
            markdown=obj.get('data',{}).get('text/markdown')
            if markdown and markdown.startswith('# PZT atomiser methods'):
                obj.data['text/markdown']=(out/'methods/METHODS.md').read_text(encoding='utf-8')
    names=['mean_v','mean_u','rms_v','rms_u','covariance','correlation','sampling_stations',
           'profile_mean_v','profile_rms_v','profile_rms_u','profile_covariance',
           'v_pdf','u_pdf','v_cdf','u_cdf','frame_879_v','frame_879_u','frame_879_quiver',
           'magnitude_mean_vector','mean_instantaneous_velocity_magnitude','velocity_magnitude_pdf']
    previews=[o for o in b.cells[3].outputs if o.output_type=='display_data' and 'image/png' in o.get('data',{})]
    assert len(previews)==len(names)
    for obj,name in zip(previews,names):
        obj.data['image/png']=base64.b64encode((out/'figures'/(name+'.png')).read_bytes()).decode('ascii')
    refreshed=len(previews)
    # Refresh the displayed configuration, summary and validation text, keeping
    # outputs traceable to current saved files and calculation cells intact.
    for obj in b.cells[1].outputs:
        if obj.output_type=='stream':obj.text=json.dumps(cfg,indent=2)+'\n'
    for obj in b.cells[4].outputs:
        if obj.output_type=='stream':
            obj.text='Measured extrema and display limits are reported separately:\n'+(out/'methods/summary.json').read_text()+'\n'+(out/'methods/validation.json').read_text()+'\n'
    b.metadata['presentation_update']='Signed contours redrawn directly from unchanged numerical arrays; matching saved-file previews and methods refreshed.'
    nb.write(b,target);shutil.copy2(target,root/'notebooks/Spatial_maps.ipynb');shutil.copy2(target,out/'methods'/target.name)
    for relative,before in baseline.items():assert digest(out/relative)==before,relative
    np.testing.assert_array_equal(old_mean,fields['mean_v'])
    np.testing.assert_allclose([old_mean.min(),old_mean.max()],[-.031319289573753975,2.2855091645694583],atol=1e-14)
    assert not mask_counts['mean_v'] and not mask_counts['sampling_stations']
    for key in r['frames']:assert not mask_counts[f'frame_{key}_v']
    report={'regenerated_contour_sets':changed,'embedded_previews_refreshed':refreshed,
            'numerical_data_files_byte_identical':list(baseline),'measurement_configuration_unchanged':True,
            'finite_values_not_masked':True,'remaining_masked_node_counts':mask_counts,
            'zero_contour_boundary_verified':True,'instantaneous_shared_limits':instant,
            'mean_axial_range_m_s':[float(old_mean.min()),float(old_mean.max())]}
    write(out/'methods/full_field_restoration_verification.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
