"""Provisional calibrated PIV: paired statistics, flag sensitivity and plots.

Internal Cartesian plot basis is [transverse V, axial U], in that order.
Every physical quantity is derived from once-converted instantaneous samples.
"""
import csv
import importlib.metadata
import json
import platform
from pathlib import Path
import numpy as np
import h5py
from .plotting import contour_map, quiver_map, shared_limits, save_figure
import matplotlib.pyplot as plt


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2), encoding='utf-8')


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def matrices(cfg):
    base_c = np.asarray(cfg['component_mapping'], float)
    base_x = np.asarray(cfg['coordinate_mapping'], float)
    common = np.asarray(cfg['transform'], float)
    for m in (base_c, base_x, common):
        if m.shape != (2,2) or not np.allclose(m @ m.T, np.eye(2)):
            raise ValueError('Mappings must be 2x2 orthogonal matrices')
    return common @ base_c, common @ base_x


def mapping_description(cfg):
    cm,_=matrices(cfg)
    formula=lambda row:' + '.join(f'({cfg["velocity_scale"]*row[i]:.9g}) {name}' for i,name in enumerate(cfg['components']))
    return f'V={formula(cm[0])} and U={formula(cm[1])}'


def coordinates(f, cfg):
    _, xm = matrices(cfg)
    horizontal, vertical = np.meshgrid(f['PARAMETERS/x'][:], f['PARAMETERS/y'][:])
    base = np.stack([horizontal,vertical]).astype(float)*cfg['coordinate_scale']
    base -= np.asarray(cfg['origin'])[:,None,None]
    return np.einsum('ab,bij->aij',xm,base)


def converted(group, cfg):
    """Scale once before calculating any statistics; do not read stored moments."""
    cm, _ = matrices(cfg)
    native = np.stack([group[k][:].astype(float) for k in cfg['components']])
    good = np.isfinite(native).all(axis=0)
    values = np.einsum('ab,bij->aij',cm,native)*cfg['velocity_scale']
    values[:,~good] = np.nan
    return values


class Moments:
    def __init__(self, shape):
        self.n = np.zeros(shape, dtype=np.int64)
        self.mean = np.zeros((2,*shape))
        self.m2 = np.zeros((2,*shape))
        self.cross = np.zeros(shape)
        self.velocity_magnitude = np.zeros(shape)

    def add(self, q, mask):
        n = self.n + mask
        delta = np.where(mask[None],q-self.mean,0.)
        self.mean += np.divide(delta,n[None],out=np.zeros_like(delta),where=n[None]>0)
        self.m2 += np.where(mask[None],delta*(q-self.mean),0.)
        self.cross += np.where(mask,delta[0]*(q[1]-self.mean[1]),0.)
        self.velocity_magnitude += np.divide(np.where(mask,np.hypot(*q)-self.velocity_magnitude,0.),n,
                                out=np.zeros_like(self.velocity_magnitude),where=n>0)
        self.n = n

    def finish(self, minimum):
        valid = self.n >= minimum
        var = np.divide(self.m2,self.n[None],out=np.full_like(self.m2,np.nan),where=valid[None])
        assert np.all(var[np.isfinite(var)] >= -1e-12)
        rms = np.sqrt(np.maximum(var,0))
        cov = np.divide(self.cross,self.n,out=np.full_like(self.cross,np.nan),where=valid)
        rho = np.divide(cov,rms[0]*rms[1],out=np.full_like(cov,np.nan),where=(rms[0]*rms[1])>0)
        finite = np.isfinite(rho)
        assert np.all(abs(rho[finite]) <= 1+1e-10)
        mean = np.where(valid[None],self.mean,np.nan)
        result = dict(mean_U=mean[1],mean_V=mean[0],
                      magnitude_mean_vector=np.hypot(*mean),
                      mean_instantaneous_velocity_magnitude=np.where(valid,self.velocity_magnitude,np.nan),
                      rms_U=rms[1],rms_V=rms[0],covariance=cov,correlation=rho)
        assert np.all(result['mean_instantaneous_velocity_magnitude'][valid]+1e-12 >= result['magnitude_mean_vector'][valid])
        return result


def calculate(path, cfg):
    with h5py.File(path,'r') as f:
        keys = sorted(f['PIV FIELDS'],key=int)
        xy = coordinates(f,cfg)  # xy[0] transverse y; xy[1] axial x
        shape = xy[0].shape
        total = Moments(shape)
        groups, flags_n = {}, {}
        frames = {}
        chosen = {keys[round(q*(len(keys)-1))] for q in cfg['representative_quantiles']}
        stations = cfg['stations']
        samples = np.full((len(keys),len(stations),2),np.nan)
        station_flags = np.zeros((len(keys),len(stations)),dtype=int)
        verify = []
        for t,key in enumerate(keys):
            g = f['PIV FIELDS'][key]
            q = converted(g,cfg)
            finite = np.isfinite(q).all(axis=0)
            assert q.shape == (2,*shape)
            total.add(q,finite)
            flags = g['flag'][:]
            for code in np.unique(flags):
                k = int(code)
                if k not in groups:
                    groups[k] = Moments(shape)
                    flags_n[k] = np.zeros(shape,dtype=np.int64)
                groups[k].add(q,finite & (flags==k))
                flags_n[k] += flags==k
            for s,station in enumerate(stations):
                i,j = station['row'],station['column']
                samples[t,s] = q[:,i,j]
                station_flags[t,s] = flags[i,j]
            if key in chosen:
                frames[key] = q.copy()
            for point in cfg['verification_vectors']:
                if str(point['frame']) != key:
                    continue
                i,j = point['row'],point['column']
                verify.append(dict(frame_key=key,row_zero_based=i,column_zero_based=j,
                    exported_x=float(f['PARAMETERS/x'][j]),exported_y=float(f['PARAMETERS/y'][i]),
                    axial_x_mm=float(xy[1,i,j]),transverse_y_mm=float(xy[0,i,j]),
                    stored_us=float(g['us'][i,j]),stored_vs=float(g['vs'][i,j]),
                    axial_U_m_s=float(q[1,i,j]),transverse_V_m_s=float(q[0,i,j]),flag_code=int(flags[i,j])))
        fields = total.finish(cfg['minimum_samples'])
        flag_fields = {k:a.finish(cfg['minimum_samples']) for k,a in groups.items()}
        cm,_ = matrices(cfg)
        expected_mean = np.einsum('ab,bij->aij',cm,np.stack([f['STATISTICS/U'][:],f['STATISTICS/V'][:]]))*cfg['velocity_scale']
        # This exporter-specific check is enabled only for the supplied us/vs population.
        stored_checks = {}
        if cfg['components'] == ['us','vs']:
            for index,key in [(0,'mean_V'),(1,'mean_U')]:
                np.testing.assert_allclose(fields[key],expected_mean[index],rtol=1e-10,atol=1e-10)
                stored_checks[key] = float(np.nanmax(abs(fields[key]-expected_mean[index])))
        checks = []
        for s,station in enumerate(stations):
            for code in [None,*sorted(groups)]:
                mask = np.isfinite(samples[:,s]).all(axis=1)
                if code is not None:
                    mask &= station_flags[:,s]==code
                a = samples[mask,s]
                if len(a)<cfg['minimum_samples']:
                    continue
                ff = fields if code is None else flag_fields[code]
                i,j = station['row'],station['column']
                expected = [a[:,1].mean(),a[:,0].mean(),a[:,1].std(),a[:,0].std(),np.mean((a[:,1]-a[:,1].mean())*(a[:,0]-a[:,0].mean()))]
                actual = [ff[k][i,j] for k in ['mean_U','mean_V','rms_U','rms_V','covariance']]
                np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
                checks.append({'station':station['id'],'population':'finite_pair' if code is None else f'flag_{code}',
                               'N':len(a),'passed':True})
        np.testing.assert_array_equal(sum(a.n for a in groups.values()),total.n)
        return dict(xy=xy,fields=fields,N=total.n,frames=frames,samples=samples,
                    station_flags=station_flags,frame_keys=keys,verification=verify,
                    flag_fields=flag_fields,flag_N={k:a.n for k,a in groups.items()},
                    flag_fraction={k:a/len(keys) for k,a in flags_n.items()},
                    metadata={'frames':len(keys),'shape':list(shape),'population':'finite-pair population',
                              'N_min':int(total.n.min()),'N_max':int(total.n.max()),
                              'independent_station_checks':checks,'converted_stored_mean_checks':stored_checks,
                              'flag_counts':{str(k):int(a.n.sum()) for k,a in groups.items()}})


VELOCITY_UNIT = r'$\mathrm{m\,s^{-1}}$'
LABELS = {'mean_U':r'$\overline{U}$ ($\mathrm{m\,s^{-1}}$)',
          'mean_V':r'$\overline{V}$ ($\mathrm{m\,s^{-1}}$)',
          'magnitude_mean_vector':r'$|\overline{\mathbf{u}}|$ ($\mathrm{m\,s^{-1}}$)',
          'mean_instantaneous_velocity_magnitude':r'$\overline{|\mathbf{u}|}$ ($\mathrm{m\,s^{-1}}$)',
          'rms_U':r'$u_{\mathrm{rms}}$ ($\mathrm{m\,s^{-1}}$)',
          'rms_V':r'$v_{\mathrm{rms}}$ ($\mathrm{m\,s^{-1}}$)',
          'covariance':r"$\overline{u'v'}$ ($\mathrm{m^2\,s^{-2}}$)",
          'correlation':r'$\rho_{uv}$'}
SIGNED = {'mean_U','mean_V','covariance','correlation'}
AXES = dict(xlabel=r'$y$ ($\mathrm{mm}$)',ylabel=r'$x$ ($\mathrm{mm}$)')


def render_maps(results,cfg,out,manifest):
    limits = {}
    for key in LABELS:
        arrays = [r['fields'][key] for r in results.values()]
        if key in ['mean_U','mean_V','rms_U','rms_V']:
            arrays += [ff[key] for r in results.values() for ff in r['flag_fields'].values()]
        limits[key] = shared_limits(arrays,key in SIGNED,key=='correlation')
    ilim = {k:shared_limits([q[i] if i<2 else np.hypot(*q) for r in results.values() for q in r['frames'].values()],i<2)
            for i,k in enumerate(['instant_V','instant_U','instant_velocity_magnitude'])}
    def contour(r,z,stem,label,lim,signed=False):
        contour_map(*r['xy'],z,label=label,limits=lim,output=stem,signed=signed,
                    levels=cfg['contour_levels'],**AXES)
        manifest.append(str(stem.relative_to(out)))
    for cid,r in results.items():
        base = out/'figures'/cid
        for key,label in LABELS.items():
            contour(r,r['fields'][key],base/key,label,limits[key],key in SIGNED)
        for key,q in r['frames'].items():
            for i,k in enumerate(['instant_V','instant_U','instant_velocity_magnitude']):
                z = q[i] if i<2 else np.hypot(*q)
                lab = [r'$V$ ($\mathrm{m\,s^{-1}}$)',r'$U$ ($\mathrm{m\,s^{-1}}$)',
                       r'$|\mathbf{u}|$ ($\mathrm{m\,s^{-1}}$)'][i]
                contour(r,z,base/f'frame_{key}_{k}',lab,ilim[k],i<2)
            stem = base/f'frame_{key}_quiver'
            quiver_map(*r['xy'],*q,output=stem,velocity_unit=VELOCITY_UNIT,reference_speed=ilim['instant_velocity_magnitude'][1],
                       stride=cfg['quiver_stride'],**AXES)
            manifest.append(str(stem.relative_to(out)))
        contour(r,r['N'],base/'finite_pair_count','Finite-pair sample count N',(0,r['metadata']['frames']))
        contour(r,r['N']/r['metadata']['frames'],base/'finite_pair_fraction','Finite-pair fraction',(0,1))
        for code,ff in r['flag_fields'].items():
            diag = out/'diagnostics'/cid/f'flag_{code}'
            contour(r,r['flag_N'][code],diag/'sample_count',f'Flag {code} finite-pair count',(0,r['metadata']['frames']))
            contour(r,r['flag_fraction'][code],diag/'fraction',f'Flag code {code} fraction',(0,1))
            for key in ['mean_U','mean_V','rms_U','rms_V']:
                contour(r,ff[key],diag/key,LABELS[key],limits[key],key in SIGNED)
    write_json(out/'results/color_limits.json',{'statistics_and_flag_diagnostics':limits,'instantaneous':ilim})


def render_profiles_and_distributions(cid,r,cfg,out,manifest):
    base = out/'figures'/cid
    xy = r['xy']
    station_rows,station_records,flag_records = [],[],[]
    for s,st in enumerate(cfg['stations']):
        i,j = st['row'],st['column']
        loc = dict(station=st['id'],row_zero_based=i,column_zero_based=j,
                   axial_x_mm=float(xy[1,i,j]),transverse_y_mm=float(xy[0,i,j]))
        station_rows.append(loc)
        for code in [None,*sorted(r['flag_fields'])]:
            finite = np.isfinite(r['samples'][:,s]).all(axis=1)
            if code is not None:
                finite &= r['station_flags'][:,s]==code
            a = r['samples'][finite,s]
            pop = 'finite_pair' if code is None else f'flag_{code}'
            record = dict(**loc,population=pop,N=len(a),fraction_of_all_frames=len(a)/len(r['frame_keys']))
            for k,column in [('U',1),('V',0)]:
                record['mean_'+k] = float(a[:,column].mean()) if len(a) else np.nan
                record['rms_'+k] = float(a[:,column].std()) if len(a)>=2 else np.nan
            flag_records.append(record)
            if code is not None:
                continue
            for k,v in [('U',a[:,1]),('V',a[:,0]),('velocity_magnitude',np.hypot(*a.T))]:
                mean,std = v.mean(),v.std()
                station_records.append(dict(**loc,quantity=k,N=len(v),finite_pair_fraction=len(v)/len(r['frame_keys']),
                    mean=float(mean),median=float(np.median(v)),rms_fluctuation=float(std),
                    p05=float(np.percentile(v,5)),p25=float(np.percentile(v,25)),p75=float(np.percentile(v,75)),
                    p95=float(np.percentile(v,95)),skewness=float(np.mean(((v-mean)/std)**3)) if std>0 else np.nan,
                    excess_kurtosis=float(np.mean(((v-mean)/std)**4)-3) if std>0 else np.nan))
    write_csv(out/'results'/f'{cid}_stations.csv',station_rows)
    write_csv(out/'results'/f'{cid}_station_statistics.csv',station_records)
    write_csv(out/'results'/f'{cid}_flag_station_comparison.csv',flag_records)
    # Each line is a native measured row. No profile interpolation or shifting.
    profile_rows = []
    for key,label in LABELS.items():
        fig,ax = plt.subplots(figsize=(6.2,4.4),layout='constrained')
        data = {}
        for n,i in enumerate(cfg['profile_rows']):
            yy,xx,z = xy[0,i],xy[1,i],r['fields'][key][i]
            ax.plot(yy,z,linestyle=['-','--',':'][n%3],label=rf'$x={xx.mean():.3f}\,\mathrm{{mm}}$')
            data[f'row_{i}_transverse_y_mm'] = yy
            data[f'row_{i}_axial_x_mm'] = xx
            data[f'row_{i}_values'] = z
            for j in range(len(yy)):
                profile_rows.append(dict(quantity=key,row_zero_based=i,column_zero_based=j,
                    axial_x_mm=float(xx[j]),transverse_y_mm=float(yy[j]),value=float(z[j]),N=int(r['N'][i,j])))
        ax.set(xlabel=AXES['xlabel'],ylabel=label)
        ax.legend(fontsize=8)
        stem = base/f'profile_{key}'
        save_figure(fig,stem,data);manifest.append(str(stem.relative_to(out)))
    write_csv(out/'results'/f'{cid}_profiles.csv',profile_rows)
    render_station_map(cid,r,cfg,out)
    manifest.append(str((base/'sampling_stations').relative_to(out)))
    # Empirical histogram PDFs: no KDE bandwidth or negative-velocity_magnitude kernel tails.
    for quantity in ['U','V','velocity_magnitude','U_normalized','V_normalized']:
        values = []
        for s,st in enumerate(cfg['stations']):
            a = r['samples'][:,s]
            a = a[np.isfinite(a).all(axis=1)]
            v = np.hypot(*a.T) if quantity=='velocity_magnitude' else a[:,1 if quantity.startswith('U') else 0]
            if quantity.endswith('normalized'):
                v=(v-v.mean())/v.std() if v.std()>0 else np.array([])
            values.append(v)
        supported=[v for v in values if len(v)]
        if not supported:
            write_json(out/'results'/f'{cid}_{quantity}_omitted.json',
                       {'reason':'No finite nonzero-variance station samples for this distribution'})
            continue
        combined = np.concatenate(supported)
        lo,hi = (0.,combined.max()) if quantity=='velocity_magnitude' else (combined.min(),combined.max())
        if hi==lo:
            lo,hi=lo-.5,hi+.5
        edges=np.linspace(lo,hi,cfg['histogram_bins']+1)
        xlabel = {'U':r'$U$ ($\mathrm{m\,s^{-1}}$)','V':r'$V$ ($\mathrm{m\,s^{-1}}$)',
                  'velocity_magnitude':r'In-plane velocity magnitude, $|\mathbf{u}|$ ($\mathrm{m\,s^{-1}}$)',
                  'U_normalized':r"$u'/u_{\mathrm{rms}}$",'V_normalized':r"$v'/v_{\mathrm{rms}}$"}[quantity]
        for kind in (['pdf'] if quantity.endswith('normalized') else ['pdf','ecdf']):
            fig,ax=plt.subplots(figsize=(6.2,4.4),layout='constrained')
            data={'bin_edges':edges}
            for n,(st,v) in enumerate(zip(cfg['stations'],values)):
                if not len(v):
                    continue
                color=f'C{n}'
                style=['-','--',':','-.','-'][n%5]
                label=f'{st["id"]}, N={len(v)}'
                density,_=np.histogram(v,bins=edges,density=True)
                data[st['id']+'_samples']=v
                data[st['id']+'_histogram_pdf']=density
                if kind=='pdf':
                    ax.stairs(density,edges,label=label,color=color,linestyle=style)
                else:
                    sx=np.sort(v); probability=np.arange(1,len(v)+1)/len(v)
                    ax.step(np.r_[sx[0],sx],np.r_[0,probability],where='post',label=label,color=color,linestyle=style)
                    data[st['id']+'_ecdf_x']=sx;data[st['id']+'_ecdf_p']=probability
            ax.set(xlabel=xlabel,ylabel='Empirical CDF' if kind=='ecdf' else ('Probability density' if quantity.endswith('normalized') else r'Probability density ($\mathrm{s\,m^{-1}}$)'))
            if quantity=='velocity_magnitude':ax.set_xlim(left=0)
            if kind=='ecdf':ax.set_ylim(0,1)
            ax.legend(fontsize=8)
            stem=base/f'{quantity}_{kind}'
            save_figure(fig,stem,data);manifest.append(str(stem.relative_to(out)))


def render_station_map(cid,r,cfg,out):
    """Show exact stations using the same axial-mean limits as all other maps."""
    xy=r['xy']
    lim=json.loads((out/'results/color_limits.json').read_text())['statistics_and_flag_diagnostics']['mean_U']
    fig,ax=plt.subplots(figsize=(6.2,5),layout='constrained')
    m=ax.contourf(*xy,np.ma.masked_invalid(r['fields']['mean_U']),levels=np.linspace(*lim,cfg['contour_levels']),
                  cmap='RdBu_r',corner_mask=False,antialiased=False)
    fig.colorbar(m,ax=ax,label=LABELS['mean_U'])
    for st in cfg['stations']:
        i,j=st['row'],st['column'];y,x=xy[:,i,j]
        ax.plot(y,x,'o',mfc='white',mec='black',ms=4)
        ax.annotate(st['id'],(y,x),xytext=(5,4),textcoords='offset points',fontsize=8)
    ax.set(**AXES,aspect='equal')
    save_figure(fig,out/'figures'/cid/'sampling_stations',
                {'transverse_y_mm':xy[0],'axial_x_mm':xy[1],'mean_U':r['fields']['mean_U'],
                 'station_rows':[s['row'] for s in cfg['stations']],
                 'station_columns':[s['column'] for s in cfg['stations']], 'color_limits':lim})


def save_results(cid,r,cfg,out):
    res=out/'results'
    np.savez_compressed(res/f'{cid}_statistics.npz',transverse_y_mm=r['xy'][0],axial_x_mm=r['xy'][1],N=r['N'],**r['fields'])
    write_csv(res/f'{cid}_spatial_statistics.csv',
        (dict(transverse_y_mm=float(r['xy'][0].ravel()[i]),axial_x_mm=float(r['xy'][1].ravel()[i]),N=int(r['N'].ravel()[i]),
              **{k:float(a.ravel()[i]) for k,a in r['fields'].items()}) for i in range(r['N'].size)))
    np.savez_compressed(res/f'{cid}_station_samples.npz',converted_V_U=r['samples'],flags=r['station_flags'],
                        frame_keys=r['frame_keys'],station_ids=[s['id'] for s in cfg['stations']])
    write_csv(res/f'{cid}_vector_comparison.csv',r['verification'])
    sheet=['# Exact-vector comparison sheet\n',
           'Configured conversion: '+mapping_description(cfg)+' (m/s). Indices below are zero-based HDF5 [row,column]; frame is the exact group key, not a timestamp. Coordinates use the exported image-reference origin, not an outlet.\n',
           '| Frame | [row,column] | Exported [x,y] | Axial x / transverse y (mm) | us / vs (assumed px/pair) | U / V (m/s) | Flag |',
           '| --- | --- | --- | --- | --- | --- | --- |']
    for p in r['verification']:
        sheet.append(f'| {p["frame_key"]} | [{p["row_zero_based"]},{p["column_zero_based"]}] | [{p["exported_x"]:g},{p["exported_y"]:g}] | {p["axial_x_mm"]:.6f} / {p["transverse_y_mm"]:.6f} | {p["stored_us"]:.8g} / {p["stored_vs"]:.8g} | {p["axial_U_m_s"]:.8g} / {p["transverse_V_m_s"]:.8g} | {p["flag_code"]} |')
    sheet += ['\nFull precision is retained in the adjacent vector_comparison.csv. Compare these exact nodes against a software numerical vector readout; if the software uses one-based array indices, add one to each row/column. Displayed coordinates may use a different reference, so check the exported coordinate fields as well.',
              '\nA magnitude factor of 0.202 should apply only for pixel displacement storage. Check vertical sign using the same vector/component readout and the software axis legend, not the overall spray direction.']
    (out/f'{cid}_VECTOR_COMPARISON.md').write_text('\n'.join(sheet),encoding='utf-8')
    write_json(res/f'{cid}_metadata.json',r['metadata'])
    sensitivity=[]
    for code,ff in r['flag_fields'].items():
        n=r['flag_N'][code]
        np.savez_compressed(res/f'{cid}_flag_{code}.npz',N=n,fraction=r['flag_fraction'][code],
                            transverse_y_mm=r['xy'][0],axial_x_mm=r['xy'][1],**ff)
        record=dict(flag_code=code,total_finite_pairs=int(n.sum()),local_N_min=int(n.min()),local_N_max=int(n.max()))
        for key in ['mean_U','mean_V','rms_U','rms_V']:
            diff=ff[key]-r['fields'][key]
            record[key+'_min']=float(np.nanmin(ff[key]));record[key+'_max']=float(np.nanmax(ff[key]))
            record[key+'_max_abs_difference_from_all_pairs']=float(np.nanmax(abs(diff)))
        sensitivity.append(record)
    write_csv(res/f'{cid}_flag_sensitivity_summary.csv',sensitivity)


METHODS = '''# Provisional calibrated results — not final validated experimental results

Spatial calibration 5.05 micrometres/pixel and pulse separation 25 microseconds
are confirmed by the user. Storage units, exported component definitions,
vertical sign and flag semantics are working assumptions, not independently
verified facts. The source arrays are us/vs because their means reproduce
exporter STATISTICS/U,V. Their upstream smoothing/replacement is undocumented.
No additional filtering, spatial smoothing, gap filling or extrapolation is
applied here.

For the current configuration, V=+0.202 us and U=+0.202 vs in m/s. The factor
is applied exactly once to instantaneous arrays, and all moments are computed
from those samples. The stored vs sign is retained provisionally; upward spray
direction and decreasing exported y do not establish its sign. No ZNMF sign
reversal is inherited. Nearby original images lack an identifiable mapping to
these export vectors and do not independently resolve signs.

Coordinate convention: x is axial/upward, y is transverse/rightward; U and V
are the corresponding velocity components. Lowercase u' and v' denote their
fluctuations about the means. The instantaneous in-plane velocity vector is
bold u=(U,V). Its magnitude is sqrt(U squared + V squared). Magnitude of the
mean in-plane velocity vector is sqrt(mean(U) squared + mean(V) squared);
mean instantaneous in-plane velocity magnitude is mean(sqrt(U squared + V
squared)). These distinct quantities are supplementary; mean_U and mean_V
are the primary velocity contours in the notebook and report selection.

Horizontal plotted y is 0.00505 times exported PARAMETERS/x;
vertical plotted axial x is 0.00505 times exported PARAMETERS/y. Origin is the
exported (0,0) image reference, not an established outlet. Both coordinate
arrays retain alignment with native rows/columns. The scaled grid spacings
are 0.3232 mm horizontally and 0.0808 mm vertically; these are vector spacing,
not interrogation size or effective resolution.

The finite-pair population consists of samples with both us and vs finite.
No flag is selected or rejected. All moments, including velocity magnitude, covariance
and component RMS, use this paired population and population denominator N.
At least two pairs are required for reported moments. Correlation is masked
at zero variance. N and finite-pair fraction are not vector-validation or
concentration metrics. Code-specific counts/means/RMS are selection-sensitivity
diagnostics; no physical acceptance or rejection meanings are assigned.

Profiles follow selected native rows without shifting or interpolation.
Stations are configured exact nodes, not pooled spatial ROIs or a physical
centreline. Locations are shown on sampling_stations and in stations.csv.
Histogram PDFs use common bin edges per quantity across stations, full sample
support and unit integral. CDFs are exact empirical staircases. Normalized
fluctuations use local N-denominator RMS; zero variance is explicitly omitted.
No KDE or negative-velocity magnitude extrapolation is used. No independence-based
confidence intervals, turbulence interpretation or periodic/stochastic
decomposition are claimed. Timing beyond the pulse separation and tracked
phase remain unverified. First/middle/last numeric field keys are illustrative,
not phase-selected or certified typical fields.

Nonfinite samples, low-count statistics and zero-variance correlations are
masked. Contours use native quadrilaterals and corner_mask=False; intrinsic
contour-boundary interpolation does not increase measurement resolution.
Common full-range limits include flag populations when comparing component
mean/RMS, and all conditions. Instantaneous scales are common across selected
frames/conditions. Every figure has PDF, SVG, 600-dpi PNG and numerical NPZ;
spatial tables, station samples and verification vectors preserve underlying
values. No plot titles are used.

To revise definitions, edit config/analysis.json and rerun the notebook or CLI.
components chooses the source pair. component_mapping maps it into [V,U];
coordinate_mapping maps exported coordinates into [transverse,axial]. transform
is a common orthogonal basis change applied to BOTH coordinates and velocities.
A later basis reflection therefore also changes arrows and recomputed signed
covariance consistently. If evidence identifies an exporter-only component
sign mismatch with already-reflected coordinates, correct component_mapping
and document why coordinate_mapping differs. Do not reflect an already-correct
coordinate twice. All statistics, samples, profiles and captions regenerate.
No derivatives are computed in this run; any future derivatives must use this
same transformed Cartesian basis. General rotations turn row profiles into
oblique cuts: the pipeline rejects such row-profile geometry rather than
mislabeling it as constant axial position.
'''


def run(config_path):
    config_path=Path(config_path).resolve();root=config_path.parents[2]
    cfg=json.loads(config_path.read_text());out=root/cfg['output']
    expected=cfg['confirmed_spatial_calibration_um_per_pixel']*1e-6/cfg['confirmed_pulse_separation_seconds']
    if not np.isclose(cfg['velocity_scale'],expected):
        if not cfg.get('verified_scale_override_reason'):
            raise ValueError('Scale differs from displacement calibration; document verified_scale_override_reason')
    for folder in ['results','figures','diagnostics','notebooks']:(out/folder).mkdir(parents=True,exist_ok=True)
    write_json(out/'configuration_used.json',cfg)
    cm,xm=matrices(cfg)
    actual=mapping_description(cfg)
    methods=METHODS.replace('V=+0.202 us and U=+0.202 vs',actual)
    methods+='\n## Configuration used\n\n'+ '\n'.join('- '+s for s in cfg['working_assumptions'])
    methods+='\n\nEffective coordinate matrix: '+str(xm.tolist())+'; effective component matrix: '+str(cm.tolist())+'.\n'
    (out/'METHODS.md').write_text(methods,encoding='utf-8')
    (out/'README.md').write_text('Provisional calibrated outputs. See METHODS.md and configuration_used.json. All assumptions remain explicit; these are not final validated measurements.\n',encoding='utf-8')
    results={};manifest=[]
    for condition in cfg['conditions']:
        cid=condition['id'];print('Calculating provisional finite-pair statistics:',cid,flush=True)
        r=calculate(root/condition['path'],cfg);results[cid]=r
        for i in cfg['profile_rows']:
            if not np.allclose(r['xy'][1,i],r['xy'][1,i,0]):
                raise ValueError('Selected row is not constant axial position after coordinate mapping')
        save_results(cid,r,cfg,out)
        print('Independent station/flag and stored-mean checks passed.',flush=True)
    render_maps(results,cfg,out,manifest)
    for cid,r in results.items():render_profiles_and_distributions(cid,r,cfg,out,manifest)
    write_json(out/'results/figure_manifest.json',manifest)
    versions={k:importlib.metadata.version(k) for k in ['numpy','h5py','matplotlib','nbformat','nbclient','ipykernel']}
    versions['python']=platform.python_version();write_json(out/'results/software_versions.json',versions)
    assumption='Provisional results; '+actual+' (m/s). '+ ' '.join(cfg['working_assumptions'])
    convention=r"Coordinate convention: $x$ is axial/upward and $y$ transverse/rightward; $U,V$ are the corresponding components, $\mathbf{u}=(U,V)$, $u'=U-\overline{U}$ and $v'=V-\overline{V}$."
    captions=['# Figure captions\n',convention+'\n',assumption+'\n']
    for stem in manifest:
        captions.append(f'**{stem}**: {caption_description(stem)}\n')
    (out/'CAPTIONS.md').write_text('\n'.join(captions),encoding='utf-8')
    write_latex_captions(out,cfg,manifest)
    write_report_selection(out,cfg,manifest)
    validate_outputs(out,results,cfg)
    write_findings(out,results,cfg)
    print(f'Exported and validated {len(manifest)} figures as PDF/SVG/600-dpi PNG plus NPZ in {out}',flush=True)
    return results


def caption_description(stem):
    name=Path(stem).name
    prefix='Native-grid transverse profile at the configured axial stations. ' if name.startswith('profile_') else ''
    key=name.removeprefix('profile_')
    descriptions={
        'mean_U':r'Mean axial velocity, $\overline{U}$.',
        'mean_V':r'Mean transverse velocity, $\overline{V}$.',
        'rms_U':r'Axial RMS fluctuation, $u_{\mathrm{rms}}=\sqrt{\overline{u\prime{}^2}}$.',
        'rms_V':r'Transverse RMS fluctuation, $v_{\mathrm{rms}}=\sqrt{\overline{v\prime{}^2}}$.',
        'covariance':r'Signed cross-covariance, $\overline{u\prime v\prime}$.',
        'correlation':r'Component correlation, $\rho_{uv}$.',
        'magnitude_mean_vector':r'Supplementary magnitude of the mean in-plane velocity vector, $|\overline{\mathbf{u}}|=\sqrt{\overline{U}^{2}+\overline{V}^{2}}$.',
        'mean_instantaneous_velocity_magnitude':r'Supplementary mean instantaneous in-plane velocity magnitude, $\overline{|\mathbf{u}|}=\overline{\sqrt{U^{2}+V^{2}}}$.',
    }
    if key in descriptions:
        result=prefix+descriptions[key]
    elif name.startswith('frame_'):
        frame=name.split('_')[1]
        quantity=(r'instantaneous $U$' if name.endswith('instant_U') else
                  r'instantaneous $V$' if name.endswith('instant_V') else
                  r'instantaneous in-plane velocity magnitude $|\mathbf{u}|$' if name.endswith('instant_velocity_magnitude') else
                  'instantaneous velocity vectors, every fourth native grid node')
        result=f'Field key {frame}: {quantity}. Common instantaneous scales are used.'
    elif '_pdf' in name or '_ecdf' in name:
        quantity=('in-plane velocity magnitude' if 'velocity_magnitude' in name else
                  'normalized component fluctuations' if 'normalized' in name else
                  r'$U$' if name.startswith('U') else r'$V$')
        result=f'Separate native-station distributions of {quantity}; histogram PDFs use full support and common bins, and CDFs are exact empirical staircases.'
        if 'velocity_magnitude' in name:result='Supplementary '+result[0].lower()+result[1:]
    elif name=='sampling_stations':
        result=r'Configured sampling nodes on the primary $\overline{U}$ contour; each PDF uses one native node.'
    else:
        result='Finite-pair sample support on the native grid.'
    if 'flag_' in str(stem):
        result+=' Flag-code population sensitivity only; code meanings are unverified and low-count statistics are masked.'
    return result


def write_report_selection(out,cfg,manifest):
    primary=[str(Path('figures')/c['id']/key) for c in cfg['conditions'] for key in cfg['primary_velocity_contours']]
    supplementary=[stem for stem in manifest if 'velocity_magnitude' in stem or 'magnitude_mean_vector' in stem]
    support=[s for s in manifest if s.startswith('figures') and s not in primary and s not in supplementary]
    selection={'primary_velocity_contours':primary,'supporting_figures':support,'supplementary_velocity_magnitude_figures':supplementary}
    assert all(s in manifest for s in primary)
    write_json(out/'results/report_figure_selection.json',selection)
    text=['# Report figure selection\n','Primary velocity contours are the two mean components. All figures retain provisional status.\n']
    for heading,stems in selection.items():
        text.append('## '+heading.replace('_',' ')+'\n')
        for stem in stems:text.append(f'- [{Path(stem).name}]({Path(stem).as_posix()}.pdf): {caption_description(stem)}')
    (out/'REPORT_FIGURE_SELECTION.md').write_text('\n'.join(text),encoding='utf-8')


def write_latex_captions(out,cfg,manifest):
    def escape(s):
        return ''.join({'_':r'\_','%':r'\%','&':r'\&','#':r'\#',
                        '{':r'\{','}':r'\}'}.get(c,c) for c in s)
    note='Provisional analysis, not final validated results. '+mapping_description(cfg)+' in m/s. Storage units and component signs are working assumptions; the finite-pair population is used without assigning validity meanings to flag codes.'
    convention=r"Coordinates: $x$ axial/upward, $y$ transverse/rightward; $\mathbf{u}=(U,V)$, $u'=U-\overline{U}$ and $v'=V-\overline{V}$."
    lines=[r'\newcommand{\PZTProvisionalNote}{'+escape(note)+' '+convention+'}',
           '% Include this shared methods note once with the caption collection.',r'\noindent\PZTProvisionalNote\par','']
    for stem in manifest:
        lines += ['% '+stem,r'\caption{'+caption_description(stem)+'}', '']
    (out/'CAPTIONS.tex').write_text('\n'.join(lines),encoding='utf-8')


def validate_outputs(out,results,cfg):
    from PIL import Image
    manifest=json.loads((out/'results/figure_manifest.json').read_text())
    for stem in manifest:
        for suffix in ['.pdf','.svg','.png','.npz']:
            p=(out/stem).with_suffix(suffix)
            assert p.exists() and p.stat().st_size>0,p
        with Image.open((out/stem).with_suffix('.png')) as im:assert abs(im.info['dpi'][0]-600)<1
        z=np.load((out/stem).with_suffix('.npz'))
        if 'values' in z:
            finite=z['values'][np.isfinite(z['values'])]
            if finite.size:
                lo,hi=z['color_limits'];assert finite.min()>=lo-1e-10 and finite.max()<=hi+1e-10
        if 'bin_edges' in z:
            for k in z.files:
                if k.endswith('_histogram_pdf'):
                    np.testing.assert_allclose(np.sum(z[k]*np.diff(z['bin_edges'])),1.,atol=1e-12)
    for cid,r in results.items():
        table=np.genfromtxt(out/'results'/f'{cid}_spatial_statistics.csv',delimiter=',',names=True)
        for key in r['fields']:np.testing.assert_allclose(table[key],r['fields'][key].ravel(),equal_nan=True)
        for key in ['U','V']:
            for st in cfg['stations']:
                i,j=st['row'],st['column']
                s=next(s for s,a in enumerate(cfg['stations']) if a['id']==st['id'])
                a=r['samples'][:,s,1 if key=='U' else 0]
                np.testing.assert_allclose(np.nanmean(a),r['fields']['mean_'+key][i,j],atol=1e-12)
    # Synthetic common-basis reflection: coordinate and velocity signs, covariance,
    # and RMS behaviour must remain mutually consistent.
    q=np.array([[[1.,2.,3.]],[[2.,1.,5.]]])
    a=Moments((1,1));b=Moments((1,1))
    for t in range(3):
        v=q[:,:,t,None];a.add(v,np.ones((1,1),bool));b.add(v*np.array([1.,-1.])[:,None,None],np.ones((1,1),bool))
    aa,bb=a.finish(2),b.finish(2)
    np.testing.assert_allclose(bb['covariance'],-aa['covariance'])
    np.testing.assert_allclose(bb['rms_U'],aa['rms_U'])
    synthetic=dict(cfg,components=['us','vs'],component_mapping=[[1,0],[0,1]],
                   coordinate_mapping=[[1,0],[0,1]],transform=[[1,0],[0,-1]],origin=[0,0])
    test_q=converted({'us':np.array([[2.]]),'vs':np.array([[3.]])},synthetic)
    test_xy=coordinates({'PARAMETERS/x':np.array([5.]),'PARAMETERS/y':np.array([7.])},synthetic)
    np.testing.assert_allclose(test_q[:,0,0],[2*cfg['velocity_scale'],-3*cfg['velocity_scale']])
    np.testing.assert_allclose(test_xy[:,0,0],[5*cfg['coordinate_scale'],-7*cfg['coordinate_scale']])
    write_json(out/'results/validation.json',dict(figure_count=len(manifest),formats_complete=True,png_dpi=600,
        histogram_integrals=True,color_ranges_cover_data=True,csv_matches_fields=True,
        station_checks=True,stored_mean_checks=True,flag_counts_partition_finite_pairs=True,
        synthetic_reflection_covariance_check=True,synthetic_coordinate_arrow_consistency=True))


def write_findings(out,results,cfg):
    text=['# Provisional numerical findings\n','All numbers use the documented working assumptions, not independently verified exporter conventions.\n']
    for cid,r in results.items():
        text += [f'## {cid}\n',f'{r["metadata"]["frames"]} fields; local finite-pair N={r["metadata"]["N_min"]} to {r["metadata"]["N_max"]}.\n',
                 '| Quantity | Spatial minimum | Spatial maximum |','| --- | ---: | ---: |']
        for k,a in r['fields'].items():text.append(f'| {k} | {np.nanmin(a):.8g} | {np.nanmax(a):.8g} |')
        i,j=np.unravel_index(np.nanargmax(r['fields']['mean_U']),r['N'].shape)
        text.append(f'\nObserved grid-node maximum mean U={r["fields"]["mean_U"][i,j]:.8g} m/s at axial x={r["xy"][1,i,j]:.8g} mm, transverse y={r["xy"][0,i,j]:.8g} mm, row={i}, column={j}, N={r["N"][i,j]}. No outlet origin or global maximum outside measured support is claimed.\n')
        text.append('Flag-specific counts and mean/RMS differences are in flag_sensitivity_summary.csv and flag_station_comparison.csv; neither group is certified physically valid.\n')
    (out/'FINDINGS.md').write_text('\n'.join(text),encoding='utf-8')
