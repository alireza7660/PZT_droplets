"""Final Cartesian PZT analysis, preserving signed samples and display masks."""
import hashlib
import json
from pathlib import Path
import platform
import importlib.metadata
import numpy as np
import h5py
from .calibrated import calculate,write_csv,write_json
from .plotting import contour_map,quiver_map,save_figure
import matplotlib.pyplot as plt

AXES={'xlabel':r'$x$ ($\mathrm{mm}$)','ylabel':r'$y$ ($\mathrm{mm}$)'}
UNIT=r'$\mathrm{m\,s^{-1}}$'
LABELS={'mean_v':r'$\overline{v}$ ($\mathrm{m\,s^{-1}}$)',
        'mean_u':r'$\overline{u}$ ($\mathrm{m\,s^{-1}}$)',
        'magnitude_mean_vector':r'$|\overline{\mathbf{w}}|$ ($\mathrm{m\,s^{-1}}$)',
        'mean_instantaneous_velocity_magnitude':r'$\overline{|\mathbf{w}|}$ ($\mathrm{m\,s^{-1}}$)',
        'rms_u':r'$u_{\mathrm{rms}}$ ($\mathrm{m\,s^{-1}}$)',
        'rms_v':r'$v_{\mathrm{rms}}$ ($\mathrm{m\,s^{-1}}$)',
        'covariance':r"$\overline{u'v'}$ ($\mathrm{m^2\,s^{-2}}$)",
        'correlation':r'$\rho_{uv}$'}
RENAME={'mean_v':'mean_U','mean_u':'mean_V','rms_v':'rms_U','rms_u':'rms_V',
        'magnitude_mean_vector':'magnitude_mean_vector',
        'mean_instantaneous_velocity_magnitude':'mean_instantaneous_velocity_magnitude',
        'covariance':'covariance','correlation':'correlation'}
PALETTE=['#1f77b4','#ff7f0e','#2ca02c','#d62728']
AXIAL_MASK_NOTE='Negative axial velocities are masked only in this contour display and the colour scale starts at zero. All signed values are retained in the numerical data and statistical calculations; no negative value is replaced by zero or an absolute value.'


def adopt(r):
    """Rename quantities only; native arrays and [horizontal,vertical] order stay unchanged."""
    return {k:r['fields'][old] for k,old in RENAME.items()}


def limits(arrays,signed=False,nonnegative=False):
    finite=np.concatenate([np.asarray(a)[np.isfinite(a)].ravel() for a in arrays])
    if nonnegative:
        positive=finite[finite>=0]
        return [0.,float(max(positive.max() if positive.size else 0.,np.finfo(float).eps))]
    lo,hi=float(finite.min()),float(finite.max())
    if lo==hi:
        delta=max(abs(lo)*.01,1e-12);lo-=delta;hi+=delta
    return [lo,hi] if signed else [0.,max(hi,np.finfo(float).eps)]


class Exporter:
    def __init__(self,out):self.out=out;self.manifest=[]
    def paths(self,name,diagnostic=False):
        stem=self.out/('diagnostics/figures' if diagnostic else 'figures')/name
        data=self.out/('diagnostics/plot_data' if diagnostic else 'data/plot_data')/(name+'.npz')
        return stem,data
    def record(self,name,caption,diagnostic=False,masked=False):
        stem,data=self.paths(name,diagnostic)
        self.manifest.append({'id':name,'category':'diagnostic' if diagnostic else ('supplementary' if 'magnitude' in name else 'primary' if name in ['mean_v','mean_u'] else 'supporting'),
          'files':{ext:str(stem.with_suffix('.'+ext).relative_to(self.out)) for ext in ['pdf','svg','png']},
          'numerical_source':str(data.relative_to(self.out)),'caption':caption+(' '+AXIAL_MASK_NOTE if masked else ''),
          'negative_axial_display_mask':masked})
    def contour(self,name,x,y,z,label,lim,caption,signed=False,masked=False,diagnostic=False,markers=None):
        stem,data=self.paths(name,diagnostic)
        contour_map(x,y,z,label=label,limits=lim,output=stem,**AXES,signed=signed,
                    mask_negative_display=masked,data_output=data,markers=markers)
        self.record(name,caption,diagnostic,masked)
    def line(self,name,fig,data,caption):
        stem,target=self.paths(name);save_figure(fig,stem,data,data_output=target);self.record(name,caption)


def render(r,fields,cfg,out,ex):
    x,y=r['xy'];nframes=r['metadata']['frames']
    mainlimits={k:limits([z],signed=k in ['mean_u','covariance','correlation'],nonnegative=k=='mean_v') for k,z in fields.items()}
    instantaneous={
        'u':limits([q[0] for q in r['frames'].values()],signed=True),
        'v':limits([q[1] for q in r['frames'].values()],nonnegative=True),
        'velocity_magnitude':limits([np.hypot(*q) for q in r['frames'].values()])}
    captions={
        'mean_v':r'Mean axial velocity $\overline{v}$.',
        'mean_u':r'Mean transverse velocity $\overline{u}$; signed measurements are shown.',
        'rms_u':r'Transverse fluctuation RMS $u_{\mathrm{rms}}$, using denominator $N$.',
        'rms_v':r'Axial fluctuation RMS $v_{\mathrm{rms}}$, using denominator $N$.',
        'covariance':r'Signed cross-covariance $\overline{u\prime v\prime}$, using denominator $N$.',
        'correlation':r'Signed correlation $\rho_{uv}$; zero-variance locations are masked.',
        'magnitude_mean_vector':r'Magnitude of the mean in-plane velocity vector, $|\overline{\mathbf{w}}|=\sqrt{\overline{u}^{2}+\overline{v}^{2}}$.',
        'mean_instantaneous_velocity_magnitude':r'Mean instantaneous in-plane velocity magnitude, $\overline{|\mathbf{w}|}=\overline{\sqrt{u^2+v^2}}$.'}
    for k,z in fields.items():
        ex.contour(k,x,y,z,LABELS[k],mainlimits[k],captions[k],signed=k in ['mean_u','covariance','correlation'],masked=k=='mean_v')
    for key,q in r['frames'].items():
        for i,component in enumerate(['u','v','velocity_magnitude']):
            z=q[i] if i<2 else np.hypot(*q)
            label=(rf'${component}$ ({UNIT})' if i<2 else rf'$|\mathbf{{w}}|$ ({UNIT})')
            caption=f'Numeric field key {key}: '+('instantaneous transverse velocity $u$.' if i==0 else 'instantaneous axial velocity $v$.' if i==1 else 'instantaneous in-plane velocity magnitude.')+' Shared limits are calculated from the three selected frames.'
            ex.contour(f'frame_{key}_{component}',x,y,z,label,instantaneous[component],caption,signed=i==0,masked=i==1)
        stem,data=ex.paths(f'frame_{key}_quiver')
        quiver_map(x,y,*q,output=stem,**AXES,velocity_unit=UNIT,
                   reference_speed=instantaneous['velocity_magnitude'][1],stride=cfg['quiver_stride'],data_output=data,
                   component_names=('u','v'))
        ex.record(f'frame_{key}_quiver',f'Instantaneous signed velocity vectors for numeric field key {key}; every fourth native node is displayed, with a common arrow scale.')
    markers=[{'x':float(x[s['row'],s['column']]),'y':float(y[s['row'],s['column']]),'id':s['id']} for s in cfg['stations']]
    ex.contour('sampling_stations',x,y,fields['mean_v'],LABELS['mean_v'],mainlimits['mean_v'],
               'Configured PDF/CDF nodes on the mean axial-velocity field; profiles use rows 287, 159 and 32.',masked=True,markers=markers)
    ex.contour('finite_pair_count',x,y,r['N'],'Finite-pair count N',[0,nframes],'Paired finite sample count, without flag selection.',diagnostic=True)
    ex.contour('finite_pair_fraction',x,y,r['N']/nframes,'Finite-pair fraction',[0,1],'Finite-pair fraction; not a physical validation measure.',diagnostic=True)
    for code,old in r['flag_fields'].items():
        ff={k:old[v] for k,v in RENAME.items()}
        ex.contour(f'flag_{code}_count',x,y,r['flag_N'][code],f'Flag {code} finite-pair count',[0,nframes],f'Code {code} count; flag meanings are not assigned.',diagnostic=True)
        ex.contour(f'flag_{code}_fraction',x,y,r['flag_fraction'][code],f'Flag {code} fraction',[0,1],f'Code {code} frequency.',diagnostic=True)
        for k in ['mean_u','mean_v','rms_u','rms_v']:
            lim=limits([ff[k]],signed=k=='mean_u',nonnegative=k=='mean_v')
            ex.contour(f'flag_{code}_{k}',x,y,ff[k],LABELS[k],lim,f'Flag {code} population sensitivity: '+captions[k]+' Diagnostic limits are independent of ensemble limits.',signed=k=='mean_u',masked=k=='mean_v',diagnostic=True)
    write_json(out/'methods/colour_limits.json',{'main_ensemble':mainlimits,'instantaneous_selected_frames':instantaneous,'main_limits_exclude_flag_diagnostics':True})
    return mainlimits


def profiles_and_distributions(r,fields,cfg,out,ex):
    x,y=r['xy'];records=[]
    for key,label in LABELS.items():
        fig,ax=plt.subplots(figsize=(6.2,4.4),layout='constrained');data={}
        for n,i in enumerate(cfg['profile_rows']):
            values=fields[key][i]
            ax.plot(x[i],values,color=PALETTE[n%4],linestyle=['-','--',':','-.'][n%4],label=rf'$y={y[i,0]:.3f}\,\mathrm{{mm}}$')
            data[f'row_{i}_x_mm']=x[i];data[f'row_{i}_y_mm']=y[i];data[f'row_{i}_values']=values
            for j,value in enumerate(values):records.append(dict(quantity=key,row=i,column=j,x_mm=float(x[i,j]),y_mm=float(y[i,j]),value=float(value),N=int(r['N'][i,j])))
        ax.set(xlabel=AXES['xlabel'],ylabel=label);ax.legend(fontsize=8)
        ex.line('profile_'+key,fig,data,'Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.')
    write_csv(out/'data/profiles.csv',records)
    stations=[];statistics=[]
    for s,st in enumerate(cfg['stations']):
        i,j=st['row'],st['column'];loc=dict(station=st['id'],row=i,column=j,x_mm=float(x[i,j]),y_mm=float(y[i,j]))
        stations.append(loc)
        for code in [None,*sorted(r['flag_fields'])]:
            good=np.isfinite(r['samples'][:,s]).all(axis=1)
            if code is not None:good &= r['station_flags'][:,s]==code
            a=r['samples'][good,s]
            for k,z in [('u',a[:,0]),('v',a[:,1]),('velocity_magnitude',np.hypot(*a.T))]:
                m=float(z.mean()) if len(z) else np.nan;sd=float(z.std()) if len(z)>=2 else np.nan
                quant=np.percentile(z,[5,25,50,75,95]) if len(z) else [np.nan]*5
                statistics.append(dict(**loc,population='finite_pair' if code is None else f'flag_{code}',quantity=k,N=len(z),fraction=len(z)/len(r['frame_keys']),mean=m,rms_fluctuation=sd,
                  p05=float(quant[0]),p25=float(quant[1]),median=float(quant[2]),p75=float(quant[3]),p95=float(quant[4])))
    write_csv(out/'data/stations.csv',stations)
    write_csv(out/'data/station_statistics.csv',[a for a in statistics if a['population']=='finite_pair'])
    write_csv(out/'diagnostics/flag_station_statistics.csv',[a for a in statistics if a['population']!='finite_pair'])
    for quantity in ['u','v','velocity_magnitude','u_normalized','v_normalized']:
        samples=[]
        for s,st in enumerate(cfg['stations']):
            a=r['samples'][:,s];a=a[np.isfinite(a).all(axis=1)]
            z=np.hypot(*a.T) if quantity=='velocity_magnitude' else a[:,0 if quantity.startswith('u') else 1]
            if quantity.endswith('normalized'):z=(z-z.mean())/z.std() if len(z)>1 and z.std()>0 else np.array([])
            samples.append(z)
        allvalues=np.concatenate([a for a in samples if len(a)])
        lo,hi=(0.,float(allvalues.max())) if quantity=='velocity_magnitude' else (float(allvalues.min()),float(allvalues.max()))
        if lo==hi:lo-=.5;hi+=.5
        bins=np.linspace(lo,hi,cfg['histogram_bins']+1)
        label=(r"$u'/u_{\mathrm{rms}}$" if quantity=='u_normalized' else r"$v'/v_{\mathrm{rms}}$" if quantity=='v_normalized' else
               r'In-plane velocity magnitude, $|\mathbf{w}|$ ('+UNIT+')' if quantity=='velocity_magnitude' else rf'${quantity}$ ({UNIT})')
        for kind in ['pdf'] if quantity.endswith('normalized') else ['pdf','cdf']:
            fig,ax=plt.subplots(figsize=(6.2,4.4),layout='constrained');data={'bin_edges':bins}
            for s,(station,z) in enumerate(zip(cfg['stations'],samples)):
                if not len(z):continue
                color=PALETTE[s%4];style=['-','--',':','-.',(0,(3,1,1,1))][s%5]
                density,_=np.histogram(z,bins=bins,density=True)
                data[station['id']+'_samples']=z;data[station['id']+'_pdf']=density
                if kind=='pdf':ax.stairs(density,bins,label=f'{station["id"]}, N={len(z)}',color=color,linestyle=style)
                else:
                    sx=np.sort(z);p=np.arange(1,len(z)+1)/len(z)
                    ax.step(np.r_[sx[0],sx],np.r_[0,p],where='post',label=f'{station["id"]}, N={len(z)}',color=color,linestyle=style)
                    data[station['id']+'_cdf_x']=sx;data[station['id']+'_cdf_probability']=p
            ax.set(xlabel=label,ylabel='Empirical CDF' if kind=='cdf' else 'Probability density' if quantity.endswith('normalized') else r'Probability density ($\mathrm{s\,m^{-1}}$)')
            if kind=='cdf':ax.set_ylim(0,1)
            if quantity=='velocity_magnitude':ax.set_xlim(left=0)
            ax.legend(fontsize=8)
            ex.line(quantity+'_'+kind,fig,data,'Separate native-node '+('empirical CDFs' if kind=='cdf' else 'full-support empirical histogram PDFs')+'; no spatial pooling. All signed component samples are retained, including negative axial velocities.')


def save_data(r,fields,cfg,out):
    x,y=r['xy']
    np.savez_compressed(out/'data/statistics.npz',x_mm=x,y_mm=y,N=r['N'],
                        variance_u=fields['rms_u']**2,variance_v=fields['rms_v']**2,**fields)
    write_csv(out/'data/spatial_statistics.csv',(dict(x_mm=float(x.ravel()[i]),y_mm=float(y.ravel()[i]),N=int(r['N'].ravel()[i]),**{k:float(a.ravel()[i]) for k,a in fields.items()}) for i in range(x.size)))
    np.savez_compressed(out/'data/station_samples.npz',u=r['samples'][:,:,0],v=r['samples'][:,:,1],flags=r['station_flags'],frame_keys=r['frame_keys'],station_ids=[s['id'] for s in cfg['stations']])
    for frame,q in r['frames'].items():np.savez_compressed(out/'data'/f'frame_{frame}.npz',x_mm=x,y_mm=y,u=q[0],v=q[1])
    for code,old in r['flag_fields'].items():np.savez_compressed(out/'diagnostics'/f'flag_{code}.npz',x_mm=x,y_mm=y,N=r['flag_N'][code],fraction=r['flag_fraction'][code],**{k:old[v] for k,v in RENAME.items()})
    verification=[]
    for p in r['verification']:
        verification.append(dict(frame_key=p['frame_key'],row=p['row_zero_based'],column=p['column_zero_based'],exported_x=p['exported_x'],exported_y=p['exported_y'],x_mm=p['transverse_y_mm'],y_mm=p['axial_x_mm'],stored_us=p['stored_us'],stored_vs=p['stored_vs'],u_m_s=p['transverse_V_m_s'],v_m_s=p['axial_U_m_s'],flag=p['flag_code']))
    write_csv(out/'data/vector_comparison.csv',verification)
    # Metadata from the shared accumulator is renamed at the external boundary.
    meta=json.loads(json.dumps(r['metadata']).replace('mean_U','mean_v').replace('mean_V','mean_u'))
    write_json(out/'methods/input_and_validation_metadata.json',meta)


def documentation(r,fields,cfg,out,ex,mainlimits):
    intro=r"Coordinates are Cartesian: $x$ is horizontal/rightward and $y$ vertical/upward; $u$ is transverse and $v$ axial. The in-plane vector is $\mathbf{w}=(u,v)$; primes denote fluctuations about the respective means."
    methods=f'''# PZT atomiser methods

{intro}

Calibration is 7.25 micrometres/pixel (0.00725 mm/pixel), replacing 5.05 for
this PZT experiment only. Pulse separation is 25 microseconds. The adopted
export convention is us/vs as horizontal/vertical pixel displacements per
pair, with retained signs: u=0.290 us and v=0.290 vs (m/s). The original HDF5
is opened read-only, converted once at the instantaneous-sample stage, then
all statistics are calculated. Naming is changed without transposing,
rotating or reflecting arrays. The exported (0,0) image reference remains
the spatial origin; it is not an established atomiser outlet centre.

The finite-pair population requires both components finite. No numeric flag
code is rejected or assigned a physical validity meaning. RMS, variance and
covariance use denominator N. Undefined correlation and fewer-than-two-pair
statistics are masked. No derivative, jet-width or similarity analysis is
added. Upstream replacement is not inferred. No additional spatial smoothing,
gap filling or extrapolation is applied.

{AXIAL_MASK_NOTE} This policy applies to mean and instantaneous axial contour
maps, the station-location map and flag-specific mean-axial contour diagnostics.
Profiles, PDFs/CDFs, quivers and all saved numerical arrays retain signed values.
Transverse velocity, covariance and correlation retain signed colour scales.
Main limits are calculated from the actual ensemble and selected instantaneous
frames; flag-specific extrema do not set ensemble colour limits. Diverging
normalization centres zero consistently without imposing symmetric endpoints.
Cells touching a display-masked node are omitted (corner_mask=False).

Magnitude of the mean vector is sqrt(mean(u)^2+mean(v)^2). Mean instantaneous
in-plane velocity magnitude is mean(sqrt(u^2+v^2)); the definitions are distinct.
Mean component contours are primary; magnitude quantities are supplementary.
Native station indices and three field keys are preserved. Profiles use the
same native rows without interpolation; PDFs are separate native-node histogram
densities with unit integral and CDFs are exact empirical staircases. No frame
independence, turbulence or absolute-uncertainty claim is made.

Grid spacings are 0.464 mm horizontally and 0.116 mm vertically; these describe
vector spacing, not interrogation-window size or effective resolution.

Regenerate from the PZT_droplets workspace root:
    .\\.venv\\Scripts\\python -m PZT_Atomiser_Analysis.analysis.pipeline --config PZT_Atomiser_Analysis/config/analysis.json
'''
    (out/'methods/METHODS.md').write_text(methods,encoding='utf-8')
    captions=['# Figure captions\n',intro+'\n']
    tex=[r'% Coordinate convention is stated once for this caption collection.',r'\noindent '+intro+r'\par','']
    for entry in ex.manifest:
        captions.append('**'+entry['id']+'**: '+entry['caption']+'\n')
        tex+=['% '+entry['id'],r'\caption{'+entry['caption']+'}','']
    (out/'methods/CAPTIONS.md').write_text('\n'.join(captions),encoding='utf-8')
    (out/'methods/CAPTIONS.tex').write_text('\n'.join(tex),encoding='utf-8')
    write_json(out/'methods/figure_manifest.json',ex.manifest)
    selection={category:[e['id'] for e in ex.manifest if e['category']==category] for category in ['primary','supporting','supplementary','diagnostic']}
    write_json(out/'methods/report_figure_selection.json',selection)
    lo,hi=float(np.nanmin(fields['mean_v'])),float(np.nanmax(fields['mean_v']))
    summary={'mean_axial_velocity_measured_range_m_s':[lo,hi],
             'mean_axial_contour_colour_limits_m_s':mainlimits['mean_v'],
             'negative_mean_axial_nodes':int(np.count_nonzero(fields['mean_v']<0)),
             'figure_count':len(ex.manifest),'publication_figure_count':sum(e['category']!='diagnostic' for e in ex.manifest),
             'diagnostic_figure_count':sum(e['category']=='diagnostic' for e in ex.manifest),
             'local_N_range':[int(r['N'].min()),int(r['N'].max())],
             'selected_frame_keys':list(r['frames'])}
    i,j=np.unravel_index(np.nanargmax(fields['mean_v']),r['N'].shape)
    summary['observed_grid_node_maximum']={'row':int(i),'column':int(j),
        'x_mm':float(r['xy'][0,i,j]),'y_mm':float(r['xy'][1,i,j]),'N':int(r['N'][i,j]),'mean_v_m_s':hi}
    write_json(out/'methods/summary.json',summary)
    (out/'README.md').write_text(f'# PZT final results\n\nMean axial measured range: {lo:.9g} to {hi:.9g} m/s. Contour display range: 0 to {mainlimits["mean_v"][1]:.9g} m/s; negatives are display-masked only.\n\n{len(ex.manifest)} figures, with data links in methods/figure_manifest.json. Methods, captions and configuration are in methods/. All raw inputs and archived results are preserved.\n',encoding='utf-8')


def validate(r,fields,cfg,out,ex):
    ratio=7.25/5.05;root=Path(__file__).resolve().parents[2]
    previous=root/'PZT_Atomiser_Analysis/provisional_calibrated/results/Batch4_Cam1_statistics.npz'
    checks={}
    with np.load(previous) as old:
        for key,oldkey in RENAME.items():
            factor=1 if key=='correlation' else ratio**2 if key=='covariance' else ratio
            np.testing.assert_allclose(fields[key],old[oldkey]*factor,rtol=1e-10,atol=1e-10,equal_nan=True)
            checks[key]=float(np.nanmax(abs(fields[key]-old[oldkey]*factor)))
        np.testing.assert_array_equal(r['N'],old['N'])
        for newkey,oldkey in [('rms_u','rms_V'),('rms_v','rms_U')]:
            np.testing.assert_allclose(fields[newkey]**2,old[oldkey]**2*ratio**2,rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(r['xy'][0],old['transverse_y_mm']*ratio,atol=1e-12)
        np.testing.assert_allclose(r['xy'][1],old['axial_x_mm']*ratio,atol=1e-12)
    from PIL import Image
    for entry in ex.manifest:
        for path in entry['files'].values():assert (out/path).stat().st_size>0
        with Image.open(out/entry['files']['png']) as im:assert abs(im.info['dpi'][0]-600)<1
        z=np.load(out/entry['numerical_source'])
        if 'values' in z:
            if entry['negative_axial_display_mask']:
                np.testing.assert_array_equal(z['display_mask'],~np.isfinite(z['values'])|(z['values']<0))
                assert z['color_limits'][0]==0
            elif 'display_mask' in z:np.testing.assert_array_equal(z['display_mask'],~np.isfinite(z['values']))
        if 'bin_edges' in z:
            for key in z.files:
                if key.endswith('_pdf'):np.testing.assert_allclose(np.sum(z[key]*np.diff(z['bin_edges'])),1.,atol=1e-12)
    with np.load(out/'data/plot_data/mean_v.npz') as p:
        np.testing.assert_array_equal(p['values'],fields['mean_v'])
    table=np.genfromtxt(out/'data/spatial_statistics.csv',delimiter=',',names=True)
    for key,z in fields.items():np.testing.assert_allclose(table[key],z.ravel(),equal_nan=True)
    write_json(out/'methods/validation.json',{'calibration_ratio':ratio,'scaling_max_absolute_errors':checks,
        'counts_and_correlations_unchanged':True,'cartesian_rename_without_transpose_or_reflection':True,
        'signed_values_preserved':True,'axial_negative_display_mask_verified':True,
        'png_600_dpi':True,'all_exports_present':True,'histogram_integrals':True,'csv_matches_npz':True})


def run(config_path):
    config_path=Path(config_path).resolve();root=config_path.parents[2];cfg=json.loads(config_path.read_text())
    expected=cfg['confirmed_spatial_calibration_um_per_pixel']*1e-6/cfg['confirmed_pulse_separation_seconds']
    np.testing.assert_allclose(cfg['velocity_scale'],expected,atol=1e-15)
    out=root/cfg['output']
    for name in ['figures','data','diagnostics','methods']:(out/name).mkdir(parents=True,exist_ok=True)
    write_json(out/'methods/configuration_used.json',cfg)
    plt.rcParams.update({'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
                         'axes.spines.top':True,'axes.spines.right':True,
                         'axes.prop_cycle':plt.cycler(color=PALETTE)})
    path=root/cfg['conditions'][0]['path'];before=(path.stat().st_size,path.stat().st_mtime_ns)
    print('Recalculating from original HDF5 at 7.25 micrometres/pixel.',flush=True)
    r=calculate(path,cfg);fields=adopt(r)
    save_data(r,fields,cfg,out)
    ex=Exporter(out);mainlimits=render(r,fields,cfg,out,ex)
    profiles_and_distributions(r,fields,cfg,out,ex)
    documentation(r,fields,cfg,out,ex,mainlimits)
    validate(r,fields,cfg,out,ex)
    assert before==(path.stat().st_size,path.stat().st_mtime_ns)
    write_json(out/'methods/software_versions.json',{'python':platform.python_version(),**{k:importlib.metadata.version(k) for k in ['numpy','h5py','matplotlib','nbformat','nbclient']}})
    print('Saved',len(ex.manifest),'figures;',json.dumps(json.loads((out/'methods/summary.json').read_text())),flush=True)
    return {'Batch4_Cam1':{'fields':fields,'x':r['xy'][0],'y':r['xy'][1],'N':r['N'],'frames':r['frames']}}
