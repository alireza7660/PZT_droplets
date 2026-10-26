"""Streaming PIV statistics and actual figure generation. Run with --config."""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import shutil
import numpy as np
import h5py
from .plotting import contour_map, quiver_map, shared_limits


def read_frame(group, cfg):
    a, b = (group[k][:].astype(float) for k in cfg['components'])
    transform = np.asarray(cfg['transform'], float)
    if not np.allclose(transform @ transform.T, np.eye(2)):
        raise ValueError('Component/coordinate transform must be orthogonal')
    u = transform[0,0]*a + transform[0,1]*b
    v = transform[1,0]*a + transform[1,1]*b
    scale = cfg['velocity_scale'] if cfg['velocity_scale'] is not None else 1.
    u *= scale
    v *= scale
    mask = np.isfinite(u) & np.isfinite(v)
    if cfg['valid_flag_values'] is not None:
        if not cfg['flag_definition_source']:
            raise ValueError('Flag selection requires documented definition source')
        mask &= np.isin(group['flag'][:], cfg['valid_flag_values'])
    return np.where(mask, u, np.nan), np.where(mask, v, np.nan)


def calculate(path, cfg):
    """Welford paired moments, denominator N. At most one full frame in memory."""
    with h5py.File(path, 'r') as f:
        keys = sorted(f['PIV FIELDS'], key=int)
        nx, ny = len(f['PARAMETERS/x']), len(f['PARAMETERS/y'])
        x, y = np.meshgrid(f['PARAMETERS/x'][:], f['PARAMETERS/y'][:])
        transform = np.asarray(cfg['transform'], float)
        scale = cfg['coordinate_scale'] if cfg['coordinate_scale'] is not None else 1.
        px, py = x*scale-cfg['origin'][0], y*scale-cfg['origin'][1]
        x, y = transform[0,0]*px+transform[0,1]*py, transform[1,0]*px+transform[1,1]*py
        n = np.zeros((ny,nx), dtype=np.int64)
        mu, mv, m2u, m2v, cross, ms = [np.zeros((ny,nx)) for _ in range(6)]
        flag_counts = {}
        chosen = {keys[round(q*(len(keys)-1))] for q in cfg['representative_quantiles']}
        frames = {}
        # Independent pointwise calculation retained for validation.
        points = [(0,0), (ny//2,nx//2), (ny-1,nx-1)]
        point_samples = [[] for _ in points]
        alternate_sum = [np.zeros((ny,nx)) for _ in range(2)]
        alternate_n = np.zeros((ny,nx), dtype=int)
        for key in keys:
            g = f['PIV FIELDS'][key]
            u, v = read_frame(g, cfg)
            assert u.shape == (ny,nx)
            good = np.isfinite(u) & np.isfinite(v)
            newn = n + good
            du, dv = np.where(good,u-mu,0.), np.where(good,v-mv,0.)
            mu += np.divide(du,newn,out=np.zeros_like(mu),where=newn>0)
            mv += np.divide(dv,newn,out=np.zeros_like(mv),where=newn>0)
            m2u += np.where(good,du*(u-mu),0.)
            m2v += np.where(good,dv*(v-mv),0.)
            cross += np.where(good,du*(v-mv),0.)
            ms += np.divide(np.where(good,np.hypot(u,v)-ms,0.),newn,
                            out=np.zeros_like(ms),where=newn>0)
            n = newn
            flags = g['flag'][:]
            for value in np.unique(flags):
                flag_counts.setdefault(int(value), np.zeros_like(n))
                flag_counts[int(value)] += flags == value
            if key in chosen:
                frames[key] = (u.copy(),v.copy())
            for samples, (i,j) in zip(point_samples,points):
                if good[i,j]:
                    samples.append((u[i,j],v[i,j]))
            if 'us' in g and 'vs' in g:
                aa, bb = g['us'][:], g['vs'][:]
                gg = np.isfinite(aa)&np.isfinite(bb)
                alternate_sum[0] += np.where(gg,aa,0)
                alternate_sum[1] += np.where(gg,bb,0)
                alternate_n += gg
        valid = n >= cfg['minimum_samples']
        varu = np.divide(m2u,n,out=np.full_like(mu,np.nan),where=valid)
        varv = np.divide(m2v,n,out=np.full_like(mu,np.nan),where=valid)
        assert np.nanmin(varu) >= -1e-12 and np.nanmin(varv) >= -1e-12
        urms, vrms = np.sqrt(np.maximum(varu,0)), np.sqrt(np.maximum(varv,0))
        cov = np.divide(cross,n,out=np.full_like(mu,np.nan),where=valid)
        rho = np.divide(cov,urms*vrms,out=np.full_like(mu,np.nan),where=(urms*vrms)>0)
        assert np.nanmax(np.abs(rho)) <= 1+1e-10
        fields = dict(mean_u=np.where(valid,mu,np.nan), mean_v=np.where(valid,mv,np.nan),
                      magnitude_mean_vector=np.where(valid,np.hypot(mu,mv),np.nan),
                      mean_instantaneous_speed=np.where(valid,ms,np.nan),
                      rms_u=urms,rms_v=vrms,covariance=cov,correlation=rho)
        assert np.all(ms[valid]+1e-10 >= np.hypot(mu,mv)[valid])
        checks = []
        for (i,j), samples in zip(points,point_samples):
            a = np.asarray(samples)
            if len(a)<2:
                continue
            expected = [*a.mean(axis=0), *a.std(axis=0), np.mean((a[:,0]-a[:,0].mean())*(a[:,1]-a[:,1].mean()))]
            actual = [mu[i,j], mv[i,j],urms[i,j],vrms[i,j],cov[i,j]]
            np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
            checks.append({'row':i,'column':j,'N':len(a),'passed':True})
        comparisons = {}
        for label, candidate in [('native_u',mu),('native_v',mv),
                                  ('alternate_us',alternate_sum[0]/np.maximum(alternate_n,1)),
                                  ('alternate_vs',alternate_sum[1]/np.maximum(alternate_n,1))]:
            stored = f['STATISTICS/'+('U' if label.endswith(('u','us')) else 'V')][:]
            comparisons[label] = float(np.nanmax(np.abs(candidate-stored)))
        metadata = {'frames':len(keys),'frame_keys':[keys[0],keys[-1]],'shape':[ny,nx],
                    'grid_spacing_native':[float(np.diff(f['PARAMETERS/x'][:])[0]),float(np.diff(f['PARAMETERS/y'][:])[0])],
                    'parameters':{k:f['PARAMETERS'][k][()].item() for k in f['PARAMETERS'] if f['PARAMETERS'][k].shape==()},
                    'independent_point_checks':checks,'stored_statistics_max_absolute_difference':comparisons,
                    'flag_counts_total':{str(k):int(a.sum()) for k,a in flag_counts.items()},
                    'local_N_range':[int(n.min()),int(n.max())],
                    'finite_pair_fraction':float(n.sum()/(len(keys)*nx*ny))}
        return dict(x=x,y=y,fields=fields,N=n,frames=frames,metadata=metadata,
                    flag_fractions={str(k):a/len(keys) for k,a in flag_counts.items()})


def plot_results(results, cfg, out):
    physical = cfg['mapping_verified']
    component_u, component_v = ('transverse V','axial U') if physical else ('component u','component v')
    vu = cfg['velocity_unit'] or 'native units; uncalibrated'
    cu = cfg['coordinate_unit'] or 'native units; uncalibrated'
    xl = f'{"Transverse position y" if physical else "Native x"} ({cu})'
    yl = f'{"Axial position x" if physical else "Native y"} ({cu})'
    labels = {'mean_u':f'Mean {component_u} ({vu})','mean_v':f'Mean {component_v} ({vu})',
              'magnitude_mean_vector':f'Magnitude of mean velocity vector ({vu})',
              'mean_instantaneous_speed':f'Mean instantaneous speed ({vu})',
              'rms_u':f'{component_u.capitalize()} RMS fluctuation ({vu})',
              'rms_v':f'{component_v.capitalize()} RMS fluctuation ({vu})',
              'covariance':f'Cross-covariance ({"("+vu+") squared"})',
              'correlation':r'Correlation $\rho_{uv}$ (dimensionless)'}
    signed = {'mean_u','mean_v','covariance','correlation'}
    limits = {k:shared_limits([r['fields'][k] for r in results.values()],k in signed,k=='correlation') for k in labels}
    # All representative fields and all conditions share the same instantaneous scales.
    instant_limits = {}
    for idx,k in enumerate(['instant_u','instant_v','instant_speed']):
        arrays = [uv[idx] if idx<2 else np.hypot(*uv) for r in results.values() for uv in r['frames'].values()]
        instant_limits[k] = shared_limits(arrays,idx<2)
    ref = instant_limits['instant_speed'][1]
    manifest = []
    for cid,r in results.items():
        base = out/'diagnostics'/'uncalibrated'/cid if not physical or not cfg['velocity_unit'] or not cfg['coordinate_unit'] else out/'figures'/cid
        args = dict(xlabel=xl,ylabel=yl)
        for k,label in labels.items():
            stem = base/k
            contour_map(r['x'],r['y'],r['fields'][k],label=label,limits=limits[k],
                        output=stem,signed=k in signed,levels=cfg['contour_levels'],**args)
            manifest.append(str(stem.relative_to(out)))
        for key,(u,v) in r['frames'].items():
            for k,z,label in [('instant_u',u,f'Instantaneous {component_u} ({vu})'),
                               ('instant_v',v,f'Instantaneous {component_v} ({vu})'),
                               ('instant_speed',np.hypot(u,v),f'Instantaneous speed ({vu})')]:
                stem = base/f'frame_{key}_{k}'
                contour_map(r['x'],r['y'],z,label=label,limits=instant_limits[k],output=stem,
                            signed=k!='instant_speed',levels=cfg['contour_levels'],**args)
                manifest.append(str(stem.relative_to(out)))
            stem = base/f'frame_{key}_quiver'
            quiver_map(r['x'],r['y'],u,v,output=stem,velocity_unit=vu,
                       reference_speed=ref,stride=cfg['quiver_stride'],**args)
            manifest.append(str(stem.relative_to(out)))
        quality = {'finite_pair_count':(r['N'],f'Finite paired sample count N (of {r["metadata"]["frames"]})',(0,r['metadata']['frames'])),
                   'finite_pair_fraction':(r['N']/r['metadata']['frames'],'Finite paired fraction (quality validity unverified)',(0,1))}
        quality.update({f'flag_{k}_fraction':(z,f'Flag code {k} fraction (meaning unverified)',(0,1)) for k,z in r['flag_fractions'].items()})
        for k,(z,label,lim) in quality.items():
            stem = base/k
            contour_map(r['x'],r['y'],z,label=label,limits=lim,output=stem,**args)
            manifest.append(str(stem.relative_to(out)))
    (out/'results/color_limits.json').write_text(json.dumps({'statistics':limits,'instantaneous':instant_limits},indent=2))
    (out/'results/figure_manifest.json').write_text(json.dumps(manifest,indent=2))
    return manifest


def run(config_path):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text())
    root = config_path.parents[2]
    out = root/cfg['output']
    for d in ['results','figures','diagnostics','Report_Files/Images']:
        (out/d).mkdir(parents=True,exist_ok=True)
    results = {}
    for condition in cfg['conditions']:
        cid = condition['id']
        print('Calculating',cid,flush=True)
        r = calculate(root/condition['path'],cfg)
        results[cid] = r
        np.savez_compressed(out/'results'/f'{cid}_statistics.npz',x=r['x'],y=r['y'],N=r['N'],**r['fields'],
                            **{'flag_'+k+'_fraction':z for k,z in r['flag_fractions'].items()})
        (out/'results'/f'{cid}_metadata.json').write_text(json.dumps(r['metadata'],indent=2))
        with (out/'results'/f'{cid}_spatial_statistics.csv').open('w',newline='') as fp:
            writer = csv.writer(fp)
            writer.writerow(['x','y','N',*r['fields']])
            writer.writerows(zip(r['x'].ravel(),r['y'].ravel(),r['N'].ravel(),*[z.ravel() for z in r['fields'].values()]))
        print('Validated',r['metadata']['independent_point_checks'],flush=True)
    manifest = plot_results(results,cfg,out)
    versions = {k:importlib.metadata.version(k) for k in ['numpy','h5py','matplotlib','nbformat','nbclient','ipykernel']}
    (out/'results/software_versions.json').write_text(json.dumps(versions,indent=2))
    print(f'Saved {len(manifest)} figures in PDF, SVG, 600-dpi PNG and NPZ.',flush=True)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',default='PZT_Atomiser_Analysis/config/analysis.json')
    run(parser.parse_args().config)
