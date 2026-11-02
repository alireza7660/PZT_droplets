"""Rigid plume-aligned measurement copies; never imports/runs interaction models."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
from .plotting import plt,contour_map

ROOT=Path(__file__).resolve().parents[2]
PROJECT=ROOT/'PZT_Atomiser_Analysis'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def json_out(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def csv_out(path,columns):
    with Path(path).open('w',newline='',encoding='utf-8') as stream:
        w=csv.writer(stream);w.writerow(columns);w.writerows(zip(*columns.values()))


def export(fig,out,stem):
    for ext in ('pdf','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)


def extract_centres(z,cfg):
    x=z['x_mm'][0];y=z['y_mm'][:,0];centres=[];dx=float(np.median(np.diff(x)))
    for i in range(len(y)):
        finite=np.isfinite(z['mean_u'][i])&np.isfinite(z['mean_v'][i])&(z['N'][i]>=cfg['minimum_finite_pair_count'])
        if not finite.any():continue
        j=int(np.argmax(np.where(finite,z['mean_v'][i],-np.inf)));peak=float(z['mean_v'][i,j])
        if j==0 or j==len(x)-1 or not finite[j-1:j+2].all() or peak<cfg['minimum_positive_peak_m_s']:continue
        a,b,c=np.polyfit(x[j-1:j+2]-x[j],z['mean_v'][i,j-1:j+2],2)
        if a>=0:continue
        delta=-b/(2*a)
        if abs(delta)>dx:continue
        left=right=j
        while left>0 and finite[left-1] and z['mean_v'][i,left-1]>=.5*peak:left-=1
        while right<len(x)-1 and finite[right+1] and z['mean_v'][i,right+1]>=.5*peak:right+=1
        nodes=right-left+1
        if nodes<cfg['minimum_connected_halfmax_native_nodes'] or left==0 or right==len(x)-1:continue
        indices=np.arange(left,right+1)
        # Secondary estimator only checks method sensitivity; never recentres rows.
        weights=z['mean_v'][i,indices]
        centroid=float(np.sum(x[indices]*weights)/np.sum(weights))
        centres.append(dict(native_row=i,peak_native_column=j,y_lab_mm=float(y[i]),x_peak_mm=float(x[j]+delta),
            x_halfmax_lobe_centroid_mm=centroid,peak_v_m_s=peak,halfmax_native_nodes=nodes,
            halfmax_node_span_mm=float(x[right]-x[left]),N_at_peak=int(z['N'][i,j])))
    return sorted(centres,key=lambda r:r['y_lab_mm'])


def fit_axis(centres,interval,key='x_peak_mm'):
    selected=[r for r in centres if interval[0]<=r['y_lab_mm']<=interval[1]]
    if len(selected)<5:raise ValueError('Insufficient supported centres in fit interval')
    y=np.array([r['y_lab_mm'] for r in selected]);x=np.array([r[key] for r in selected])
    slope,intercept=np.polyfit(y,x,1);res=x-(intercept+slope*y)
    return dict(requested_lower_y_mm=interval[0],requested_upper_y_mm=interval[1],
        actual_lower_y_mm=float(y.min()),actual_upper_y_mm=float(y.max()),centre_count=len(y),
        slope_dx_dy=float(slope),intercept_x_at_lab_y0_mm=float(intercept),
        angle_right_of_lab_vertical_deg=float(np.degrees(np.arctan(slope))),
        residual_RMS_x_mm=float(np.sqrt(np.mean(res*res))),residual_max_abs_x_mm=float(np.max(abs(res))),
        halfmax_nodes_min=min(r['halfmax_native_nodes'] for r in selected),halfmax_nodes_max=max(r['halfmax_native_nodes'] for r in selected),estimator=key)


def sample_native(z,points_lab):
    """Two-component bilinear samples; invalid/outside cells remain NaN."""
    x=z['x_mm'][0];order=np.argsort(z['y_mm'][:,0]);y=z['y_mm'][order,0]
    u=z['mean_u'][order];v=z['mean_v'][order];points=np.asarray(points_lab)
    answer=np.full(points.shape,np.nan)
    for k,(xx,yy) in enumerate(points):
        if not (x[0]<=xx<=x[-1] and y[0]<=yy<=y[-1]):continue
        j=min(max(int(np.searchsorted(x,xx,side='right'))-1,0),len(x)-2)
        i=min(max(int(np.searchsorted(y,yy,side='right'))-1,0),len(y)-2)
        block_u=u[i:i+2,j:j+2];block_v=v[i:i+2,j:j+2]
        if not (np.isfinite(block_u).all() and np.isfinite(block_v).all()):continue
        wx=(xx-x[j])/(x[j+1]-x[j]);wy=(yy-y[i])/(y[i+1]-y[i])
        w=np.array([[(1-wx)*(1-wy),wx*(1-wy)],[(1-wx)*wy,wx*wy]])
        answer[k]=[np.sum(w*block_u),np.sum(w*block_v)]
    return answer


def tensor_rotate(z,Q):
    C=np.empty(z['mean_u'].shape+(2,2));C[...,0,0]=z['variance_u'];C[...,1,1]=z['variance_v']
    C[...,0,1]=C[...,1,0]=z['covariance']
    transformed=np.einsum('ab,...bc,dc->...ad',Q,C,Q)
    return C,transformed


def half_widths(n,v):
    valid=np.isfinite(v)
    if not valid.any():return dict(peak_n_mm=None,left_half_width_mm=None,right_half_width_mm=None,width_asymmetry=None)
    k=int(np.nanargmax(v));peak=v[k]
    crossings=[]
    for step in (-1,1):
        last=k;found=None
        for j in range(k+step,len(n) if step>0 else -1,step):
            if not np.isfinite(v[j]):break
            if v[last]>=.5*peak and v[j]<=.5*peak:
                found=n[last]+(.5*peak-v[last])*(n[j]-n[last])/(v[j]-v[last]);break
            last=j
        crossings.append(found)
    if None in crossings:return dict(peak_n_mm=float(n[k]),left_half_width_mm=None,right_half_width_mm=None,width_asymmetry=None)
    left=n[k]-crossings[0];right=crossings[1]-n[k]
    return dict(peak_n_mm=float(n[k]),left_half_width_mm=float(left),right_half_width_mm=float(right),width_asymmetry=float((right-left)/(right+left)))


def run(config_path=None):
    cfg_path=Path(config_path or PROJECT/'config/plume_alignment.json');cfg=json.loads(cfg_path.read_text())
    source=ROOT/cfg['statistics'];out=ROOT/cfg['output']
    for folder in ('data','figures','reference','validation'):(out/folder).mkdir(parents=True,exist_ok=True)
    final=source.parent.parent
    protected={str(p.resolve()):sha(p) for p in final.rglob('*') if p.is_file() and out not in p.resolve().parents and p.suffix.lower() in ('.npz','.csv','.pdf','.png','.svg')}
    with np.load(source) as a:z={k:a[k].copy() for k in a.files}
    source_cfg=json.loads((ROOT/cfg['source_configuration']).read_text())
    np.testing.assert_allclose([source_cfg['coordinate_scale'],source_cfg['velocity_scale']],[.00535,.214],rtol=0,atol=1e-15)
    centres=extract_centres(z,cfg);fit=fit_axis(centres,cfg['fit_interval_lab_y_mm'])
    sensitivities=[fit_axis(centres,interval) for interval in cfg['sensitivity_intervals_lab_y_mm']]
    method_sensitivity=fit_axis(centres,cfg['fit_interval_lab_y_mm'],'x_halfmax_lobe_centroid_mm')
    b=fit['slope_dx_dy'];angle=np.arctan(b);c=np.cos(angle);s=np.sin(angle);Q=np.array([[c,-s],[s,c]])
    y_ref=fit['actual_lower_y_mm'];origin=np.array([fit['intercept_x_at_lab_y0_mm']+b*y_ref,y_ref])
    positions=np.stack([z['x_mm'],z['y_mm']],axis=-1);vel=np.stack([z['mean_u'],z['mean_v']],axis=-1)
    aligned=(positions-origin)@Q.T;velocity=vel@Q.T;C,Cp=tensor_rotate(z,Q)
    if np.nanmin(np.diagonal(Cp,axis1=-2,axis2=-1))<0:raise ValueError('Negative rotated variance: inspect source covariance')
    N=z['N'].copy();n=aligned[...,0];sp=aligned[...,1];up=velocity[...,0];vp=velocity[...,1]
    arrays=dict(n_p_mm=n,s_p_mm=sp,mean_u_p=up,mean_v_p=vp,N=N,
        variance_u_p=Cp[...,0,0],variance_v_p=Cp[...,1,1],covariance_u_p_v_p=Cp[...,0,1],
        rms_u_p=np.sqrt(Cp[...,0,0]),rms_v_p=np.sqrt(Cp[...,1,1]),
        correlation_u_p_v_p=Cp[...,0,1]/np.sqrt(Cp[...,0,0]*Cp[...,1,1]),
        magnitude_mean_vector=z['magnitude_mean_vector'],mean_instantaneous_velocity_magnitude=z['mean_instantaneous_velocity_magnitude'],
        original_lab_x_mm=z['x_mm'],original_lab_y_mm=z['y_mm'],original_finite_pair_mask=np.isfinite(vel).all(axis=-1))
    np.savez_compressed(out/'data/statistics_plume_aligned.npz',**arrays)
    rows,cols=np.indices(up.shape)
    csv_out(out/'data/plume_aligned_mean_fields.csv',dict(native_row=rows.ravel(),native_column=cols.ravel(),
        original_lab_x_mm=z['x_mm'].ravel(),original_lab_y_mm=z['y_mm'].ravel(),n_p_mm=n.ravel(),s_p_mm=sp.ravel(),
        mean_u_p_m_s=up.ravel(),mean_v_p_m_s=vp.ravel(),N=N.ravel()))
    for row in centres:
        row['in_primary_fit']=cfg['fit_interval_lab_y_mm'][0]<=row['y_lab_mm']<=cfg['fit_interval_lab_y_mm'][1]
        row['fitted_x_mm']=fit['intercept_x_at_lab_y0_mm']+b*row['y_lab_mm'];row['residual_x_mm']=row['x_peak_mm']-row['fitted_x_mm']
        transformed=(np.array([row['x_peak_mm'],row['y_lab_mm']])-origin)@Q.T
        row['centre_n_p_mm']=float(transformed[0]);row['centre_s_p_mm']=float(transformed[1])
    csv_out(out/'data/extracted_centres.csv',{k:[r[k] for r in centres] for k in centres[0]})
    csv_out(out/'data/fit_interval_sensitivity.csv',{k:[r[k] for r in sensitivities] for k in sensitivities[0]})
    # Fixed fitted-axis samples: native lab y planes, different x from the original fixed column.
    y_axis=np.sort(z['y_mm'][:,0]);lab_axis=np.c_[fit['intercept_x_at_lab_y0_mm']+b*y_axis,y_axis]
    axis_velocity=sample_native(z,lab_axis)@Q.T;axis_points=(lab_axis-origin)@Q.T
    axis_data=dict(s_p_mm=axis_points[:,1],n_p_mm=axis_points[:,0],lab_x_mm=lab_axis[:,0],lab_y_mm=lab_axis[:,1],
                   mean_u_p_m_s=axis_velocity[:,0],mean_v_p_m_s=axis_velocity[:,1],finite_pair=np.isfinite(axis_velocity).all(axis=1))
    csv_out(out/'data/fitted_axis_axial_profile.csv',axis_data);np.savez_compressed(out/'data/fitted_axis_axial_profile.npz',**axis_data)
    dx=float(np.median(np.diff(z['x_mm'][0])));stations=[];section_arrays=[]
    nn=np.arange(np.floor(n.min()/dx),np.ceil(n.max()/dx)+1)*dx
    for station in cfg['transverse_stations_plume_s_mm']:
        points=np.c_[nn,np.full(len(nn),station)]@Q+origin
        vv=sample_native(z,points)@Q.T;metrics=half_widths(nn,vv[:,1]);stations.append(dict(s_p_mm=station,
            axis_lab_x_mm=float((np.array([0,station])@Q+origin)[0]),axis_lab_y_mm=float((np.array([0,station])@Q+origin)[1]),**metrics))
        section_arrays.append(dict(n_p_mm=nn,s_p_mm=np.full(len(nn),station),lab_x_mm=points[:,0],lab_y_mm=points[:,1],
                                   mean_u_p_m_s=vv[:,0],mean_v_p_m_s=vv[:,1],finite_pair=np.isfinite(vv).all(axis=1)))
    csv_out(out/'data/transverse_profiles.csv',{k:np.concatenate([a[k] for a in section_arrays]) for k in section_arrays[0]})
    csv_out(out/'data/transverse_profile_metrics.csv',{k:[r[k] for r in stations] for k in stations[0]})
    # Preserve original fixed-column exports byte-for-byte as separate references.
    for ext in ('csv','npz'):shutil.copy2(final/'data'/('mean_v_vertical_profile.'+ext),out/'reference'/('original_fixed_column_profile.'+ext))
    with np.load(final/'data/mean_v_vertical_profile.npz') as a:original_profile={k:a[k].copy() for k in a.files}
    profile_rows=original_profile['native_row'];profile_cols=original_profile['native_column']
    fixed_s=sp[profile_rows,profile_cols];fixed_vp=vp[profile_rows,profile_cols]
    csv_out(out/'data/original_fixed_column_in_aligned_coordinates.csv',dict(native_row=profile_rows,native_column=profile_cols,
        lab_x_mm=z['x_mm'][profile_rows,profile_cols],lab_y_mm=z['y_mm'][profile_rows,profile_cols],
        n_p_mm=n[profile_rows,profile_cols],s_p_mm=fixed_s,original_mean_lab_v_m_s=z['mean_v'][profile_rows,profile_cols],mean_v_p_m_s=fixed_vp))
    gravity=Q@np.asarray(cfg['gravity_lab_m_s2'])
    metadata=dict(primary_fit=fit,fit_sensitivity=sensitivities,secondary_centre_estimator_fit=method_sensitivity,
        orthonormal_matrix_lab_to_plume=Q.tolist(),inverse_matrix_plume_to_lab=Q.T.tolist(),reference_point_lab_mm=origin.tolist(),
        reference_description='Fitted-axis intersection with lowest accepted native y in primary interval; not atomiser surface',
        formulas='[n_p,s_p]=Q*([x_lab,y_lab]-reference); [u_p,v_p]=Q*[u_lab,v_lab]; C_p=Q*C_lab*Q^T',
        gravity_lab_m_s2=cfg['gravity_lab_m_s2'],gravity_plume_m_s2=gravity.tolist(),axial_label='s_p: plume-aligned axial, not laboratory vertical',
        coordinate_units='mm',velocity_units='m s^-1',variance_covariance_units='m^2 s^-2',N_meaning='Original finite paired image-field sample count, not droplet concentration',
        frame_interpretation='Coordinate representation of the measured flow only; no claim of identical flow after physical remounting',
        coordinates_resampled_for_contours=False,profile_interpolation='Four-finite-corner bilinear samples on fixed straight axis and fixed s_p cuts, no extrapolation',
        original_fixed_column_profile_preserved=True,image_evidence='Nearby EVOS microscope TIFFs and assembly illustration do not establish orientation for this PIV acquisition; see validation/image_inventory.json')
    json_out(out/'transformation.json',metadata);json_out(out/'configuration_used.json',cfg)
    manifest=render(out,z,arrays,centres,fit,origin,Q,axis_data,section_arrays,original_profile,fixed_s,fixed_vp)
    json_out(out/'report_manifest.json',dict(figures=manifest,text=['README.md','methods.tex','captions.tex'],interaction_predictions_excluded=True))
    if np.nanmax(abs(axis_points[:,0]))>1e-12:raise AssertionError('Fitted-axis samples are not at transverse zero')
    checks=verify(z,arrays,Q,origin,C,Cp,positions,velocity,gravity,protected)
    checks['protected_files_hash_verified']=len(protected);json_out(out/'validation/numerical_verification.json',checks)
    json_out(out/'validation/protected_file_hashes.json',protected)
    json_out(out/'source_manifest.json',dict(statistics_path=str(source.resolve()),statistics_sha256=sha(source),
        processing_configuration=str((ROOT/cfg['source_configuration']).resolve()),processing_configuration_sha256=sha(ROOT/cfg['source_configuration']),
        physical_calibration='5.35 um/pixel, 25 us, 0.214 m/s per pixel/pair already applied; no recalibration',
        component_source='Saved mean_u horizontal/right, mean_v laboratory vertical/up; rotation uses both components',
        sample_population='Original finite-pair us/vs population unchanged; N copied, no flag reclassification',
        covariance_source='Saved variance_u, variance_v and covariance, population denominator N',
        original_image_inspection='RAW_data has HDF5 only. Nearby Atomiser images/VMARF TIFF metadata identifies EVOS microscope and 0.619 um/pixel, not the linked PIV 5.35 um/pixel image sequence. Assembly.png is an illustration, not acquisition-angle evidence.'))
    from .plume_alignment_report import write_report
    write_report(out,cfg,metadata,stations,manifest,checks)
    return dict(angle_right_of_lab_vertical_deg=fit['angle_right_of_lab_vertical_deg'],fit_interval_actual_lab_y_mm=[fit['actual_lower_y_mm'],fit['actual_upper_y_mm']],
        centre_count=fit['centre_count'],fit_residual_RMS_mm=fit['residual_RMS_x_mm'],fit_residual_max_mm=fit['residual_max_abs_x_mm'],
        sensitivity_angle_range_deg=[min(r['angle_right_of_lab_vertical_deg'] for r in sensitivities),max(r['angle_right_of_lab_vertical_deg'] for r in sensitivities)],
        reference_point_lab_mm=origin.tolist(),gravity_plume_m_s2=gravity.tolist(),transverse_metrics=stations,
        protected_files_unchanged=len(protected),interaction_predictions_run=False)


def verify(z,aligned,Q,origin,C,Cp,positions,velocity,gravity,protected):
    np.testing.assert_allclose(Q@Q.T,np.eye(2),atol=1e-15);np.testing.assert_allclose(np.linalg.det(Q),1,atol=1e-15)
    restored=np.stack([aligned['n_p_mm'],aligned['s_p_mm']],axis=-1)@Q+origin
    np.testing.assert_allclose(restored,positions,rtol=1e-14,atol=1e-13)
    restored_velocity=velocity@Q
    original=np.stack([z['mean_u'],z['mean_v']],axis=-1)
    np.testing.assert_allclose(restored_velocity,original,rtol=1e-12,atol=1e-15)
    np.testing.assert_allclose(np.linalg.norm(velocity,axis=-1),np.linalg.norm(original,axis=-1),rtol=1e-14,atol=1e-15)
    restored_C=np.einsum('ba,...bc,cd->...ad',Q,Cp,Q)
    np.testing.assert_allclose(restored_C,C,rtol=1e-11,atol=1e-15)
    np.testing.assert_allclose(np.trace(Cp,axis1=-2,axis2=-1),np.trace(C,axis1=-2,axis2=-1),atol=1e-15)
    np.testing.assert_allclose(np.linalg.eigvalsh(Cp),np.linalg.eigvalsh(C),rtol=1e-10,atol=1e-15)
    np.testing.assert_array_equal(aligned['N'],z['N']);np.testing.assert_array_equal(np.isfinite(velocity).all(axis=-1),np.isfinite(original).all(axis=-1))
    np.testing.assert_allclose(Q.T@gravity,[0,-9.80665],atol=1e-14)
    # Deliberately outside points must not receive extrapolated samples.
    assert np.isnan(sample_native(z,[[z['x_mm'].min()-1,z['y_mm'].min()-1]])).all()
    for path,value in protected.items():assert sha(path)==value,path
    return dict(orthonormal_and_det_plus_one=True,position_inverse_verified=True,velocity_inverse_verified=True,
        maximum_velocity_magnitude_difference_m_s=float(np.max(abs(np.linalg.norm(velocity,axis=-1)-np.linalg.norm(original,axis=-1)))),
        covariance_tensor_inverse_and_eigenvalues_verified=True,RMS_derived_from_rotated_covariance=True,
        original_population_and_masks_preserved=True,no_native_grid_resampling=True,outside_profile_samples_remain_missing=True,
        transformed_gravity_verified=True,interaction_predictions_run=False)


def render(out,z,a,centres,fit,origin,Q,axis_data,sections,original_profile,fixed_s,fixed_vp):
    from matplotlib.colors import Normalize
    figs=[];extent=max(float(np.nanmax(abs(z['mean_v']))),float(np.nanmax(abs(a['mean_v_p']))));norm=Normalize(-extent,extent)
    levels=np.linspace(-extent,extent,41)
    def native(ax,x,y,v):
        m=ax.contourf(x,y,np.ma.masked_invalid(v),levels=levels,cmap='RdBu_r',norm=norm,corner_mask=False,antialiased=False)
        ax.set_aspect('equal');return m
    fig,ax=plt.subplots(figsize=(6.4,5.3),layout='constrained');m=native(ax,z['x_mm'],z['y_mm'],z['mean_v'])
    for flag,color,label in [(False,'0.6','Centres outside fit'),(True,'black','Centres used in fit')]:
        group=[r for r in centres if r['in_primary_fit']==flag];ax.plot([r['x_peak_mm'] for r in group],[r['y_lab_mm'] for r in group],'.',ms=2.5,color=color,label=label)
    y=np.linspace(z['y_mm'].min(),z['y_mm'].max(),200);ax.plot(fit['intercept_x_at_lab_y0_mm']+fit['slope_dx_dy']*y,y,'k--',lw=1,label='One fitted straight axis')
    ax.plot(origin[0],origin[1],'o',mfc='white',mec='black',ms=4);ax.set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$');ax.legend(frameon=False,fontsize=8)
    fig.colorbar(m,ax=ax,label=r'$\overline{v}\ (\mathrm{m\,s^{-1}})$');stem='original_field_axis_fit';export(fig,out,stem);figs.append(stem)
    for field,label,stem in [('mean_v_p',r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$','plume_aligned_mean_axial_velocity'),
                             ('mean_u_p',r'$\overline{u}_p\ (\mathrm{m\,s^{-1}})$','plume_aligned_mean_transverse_velocity')]:
        limit=extent if field=='mean_v_p' else float(np.nanmax(abs(a[field])))
        contour_map(a['n_p_mm'],a['s_p_mm'],a[field],label=label,limits=(-limit,limit),signed=True,
            output=out/'figures'/stem,xlabel=r'$n_p\ (\mathrm{mm})$',ylabel=r'$s_p\ (\mathrm{mm})$',data_output=out/'data'/(stem+'.npz'))
        figs.append(stem)
    fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
    ax.plot(axis_data['s_p_mm'],axis_data['mean_v_p_m_s'],color='#d62728',label=r'Fixed fitted axis: $n_p=0$')
    ax.set(xlabel=r'$s_p\ (\mathrm{mm})$',ylabel=r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$');ax.legend(frameon=False)
    stem='fitted_axis_axial_velocity_profile';export(fig,out,stem);figs.append(stem)
    for component,label,stem in [('mean_v_p_m_s',r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$','transverse_axial_velocity_profiles'),
                                  ('mean_u_p_m_s',r'$\overline{u}_p\ (\mathrm{m\,s^{-1}})$','transverse_transverse_velocity_profiles')]:
        fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
        for section in sections:
            core=np.isfinite(section[component])&(abs(section['n_p_mm'])<=7)
            values=np.where(core,section[component],np.nan)
            ax.plot(section['n_p_mm'],values,label=rf'$s_p={section["s_p_mm"][0]:g}$ mm')
        ax.axvline(0,color='0.5',ls=':',lw=.7);ax.set(xlabel=r'$n_p\ (\mathrm{mm})$',ylabel=label,xlim=(-7,7));ax.legend(frameon=False,fontsize=8)
        export(fig,out,stem);figs.append(stem)
    fig,axes=plt.subplots(1,2,figsize=(10,5.2),layout='constrained')
    native(axes[0],z['x_mm'],z['y_mm'],z['mean_v']);m=native(axes[1],a['n_p_mm'],a['s_p_mm'],a['mean_v_p'])
    axes[0].set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$');axes[1].set(xlabel=r'$n_p\ (\mathrm{mm})$',ylabel=r'$s_p\ (\mathrm{mm})$')
    fig.colorbar(m,ax=list(axes),label=r'$\overline{v},\ \overline{v}_p\ (\mathrm{m\,s^{-1}})$',shrink=.8)
    stem='original_vs_plume_aligned';export(fig,out,stem);figs.append(stem)
    fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
    ax.plot(fixed_s,original_profile['mean_v'],color='0.5',ls=':',label='Original fixed column, laboratory component')
    ax.plot(fixed_s,fixed_vp,'k--',label='Same fixed column, aligned component')
    ax.plot(axis_data['s_p_mm'],axis_data['mean_v_p_m_s'],color='#d62728',label='Fitted-axis line, aligned component')
    ax.set(xlabel=r'$s_p\ (\mathrm{mm})$',ylabel=r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8)
    stem='different_profile_lines_comparison';export(fig,out,stem);figs.append(stem)
    fig,ax=plt.subplots(figsize=(6.4,4.3),layout='constrained')
    ax.plot([r['y_lab_mm'] for r in centres],[r['residual_x_mm'] for r in centres],'.-',ms=2,lw=.7)
    ax.axhline(0,color='black',lw=.7);ax.axvspan(*[fit['actual_lower_y_mm'],fit['actual_upper_y_mm']],alpha=.1,color='tab:blue')
    ax.set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$x_c-x_{\mathrm{fit}}\ (\mathrm{mm})$')
    stem='straight_axis_residuals';export(fig,out,stem);figs.append(stem)
    return figs


def notebook_cells():
    import nbformat as nb
    text='''## Plume-aligned measured field — interaction prediction paused

One straight hydrodynamic axis is fitted over documented laboratory heights,
with interval sensitivity and residuals reported. A single orthonormal transform
acts on positions, velocities, covariance and gravity. The rotated native grid
is plotted without resampling. The axial coordinate s_p is plume-aligned, not
laboratory vertical; its zero is an analysis reference, not the atomiser surface.
This does not simulate physical remounting or remove genuine asymmetry.
'''
    code='''from PZT_Atomiser_Analysis.analysis.plume_alignment import run as align_plume
alignment=align_plume()
print(json.dumps(alignment,indent=2))
aligned_output=output/'Plume_Aligned'
for stem in json.loads((aligned_output/'report_manifest.json').read_text())['figures']:
    display(Image(filename=str(aligned_output/'figures'/(stem+'.png')),width=650))
'''
    cells=[nb.v4.new_markdown_cell(text),nb.v4.new_code_cell(code)]
    for c in cells:c.metadata['plume_alignment']=True
    return cells


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path)
    print(json.dumps(run(p.parse_args().config),indent=2))
