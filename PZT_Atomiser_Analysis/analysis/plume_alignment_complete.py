"""Complete native-grid plume-aligned exports and standalone report figures."""
from pathlib import Path
import json
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from .plotting import plt


def render_complete(out,z,a,centres,fit,origin,Q,axis,sections,original,fixed_y,fixed_v,cfg,metadata):
    from .plume_alignment import csv_out,json_out,sha,sample_native
    x,y=a['x_p_mm'],a['y_p_mm'];manifest=[]
    # Width contains observed half-maximum lobes, not an imposed symmetric field.
    radius=0.
    for row in centres:
        i=row['native_row'];j=row['peak_native_column'];v=z['mean_v'][i];valid=np.isfinite(v)
        left=right=j
        while left>0 and valid[left-1] and v[left-1]>=.5*v[j]:left-=1
        while right<len(v)-1 and valid[right+1] and v[right+1]>=.5*v[j]:right+=1
        radius=max(radius,float(np.max(abs(x[i,left:right+1]))))
    policy=cfg['report_crop'];step=policy['half_width_rounding_mm'];w=np.ceil((radius+policy['halfmax_margin_mm'])/step)*step
    c,s=Q[0,0],Q[1,0];xmin,xmax=float(z['x_mm'].min()),float(z['x_mm'].max());ymin,ymax=float(z['y_mm'].min()),float(z['y_mm'].max())
    # All four corners satisfy rectangular lab support. x_lab=x0+c*xp+s*yp.
    low=(ymin-origin[1]+abs(s)*w)/c;high=(ymax-origin[1]-abs(s)*w)/c
    if s>0:low=max(low,(xmin-origin[0]+c*w)/s);high=min(high,(xmax-origin[0]-c*w)/s)
    elif s<0:low=max(low,(xmax-origin[0]-c*w)/s);high=min(high,(xmin-origin[0]+c*w)/s)
    assert high>low
    corners=np.array([[-w,low],[w,low],[-w,high],[w,high]])
    lab=corners@Q+origin
    assert np.all(lab[:,0]>=xmin-1e-10) and np.all(lab[:,0]<=xmax+1e-10) and np.all(lab[:,1]>=ymin-1e-10) and np.all(lab[:,1]<=ymax+1e-10)
    native_inside=(abs(x)<=w)&(y>=low)&(y<=high)
    full_bounds=[-float(np.max(abs(x))),float(np.max(abs(x))),float(y.min()),float(y.max())]
    crop=dict(x_p_bounds_mm=[-float(w),float(w)],y_p_bounds_mm=[float(low),float(high)],lab_corners_mm=lab.tolist(),
        method=policy['method'],extracted_halfmax_max_radius_mm=radius,internal_masks_preserved=True,
        full_support_y_p_bounds_mm=[float(y.min()),float(y.max())],retained_axial_fraction=float((high-low)/np.ptp(y)),retained_fitted_axis_fraction=float((high-low)/((ymax-ymin)/c)),
        excluded_axial_mm=float(np.ptp(y)-(high-low)),native_node_fraction_inside_crop=float(native_inside.mean()),
        note='Report axes clip native contour polygons; native arrays are never cropped or resampled. Rectangle corners inside lab domain do not fill internal invalid masks.')
    metadata['report_crop']=crop;metadata['inverse_position_formula']='[x_lab,y_lab]=Q^T*[x_p,y_p]+reference';metadata['inverse_velocity_formula']='[u_lab,v_lab]=Q^T*[u_p,v_p]'
    metadata['station_y_p_mm']=cfg['transverse_stations_plume_y_mm'];json_out(out/'transformation.json',metadata)
    np.savez_compressed(out/'data/report_crop_support.npz',native_node_inside_report_rectangle=native_inside,rectangle_corners_plume_mm=corners,rectangle_corners_lab_mm=lab)
    labels={'mean_v_p':r'$\overline{v}_p\ (\mathrm{m\,s^{-1}})$','mean_u_p':r'$\overline{u}_p\ (\mathrm{m\,s^{-1}})$',
        'rms_v_p':r'$v_{p,\mathrm{rms}}\ (\mathrm{m\,s^{-1}})$','rms_u_p':r'$u_{p,\mathrm{rms}}\ (\mathrm{m\,s^{-1}})$',
        'covariance_u_p_v_p':r'$\overline{u_p^{\prime}v_p^{\prime}}\ (\mathrm{m^2\,s^{-2}})$','correlation_u_p_v_p':r'$\rho_{u_pv_p}$',
        'N':r'$N\ (\mathrm{paired\ samples})$'}
    stems={'mean_v_p':'mean_axial_velocity','mean_u_p':'mean_transverse_velocity','rms_v_p':'axial_RMS','rms_u_p':'transverse_RMS','covariance_u_p_v_p':'axial_transverse_covariance','correlation_u_p_v_p':'correlation','N':'sample_count'}
    def save(fig,stem,category,caption,data=None):
        folder=out/('figures/report' if category=='report' else 'diagnostics/figures' if category=='diagnostic' else 'diagnostics/full_support');folder.mkdir(parents=True,exist_ok=True)
        for ext in ('pdf','svg','png'):fig.savefig(folder/(stem+'.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig)
        manifest.append(dict(stem=stem,path=str((folder/(stem+'.pdf')).relative_to(out)).replace('\\','/'),category=category,caption=caption,formats=['pdf','svg','png'],PNG_dpi=600))
        if data is not None:np.savez_compressed(out/'data'/(category+'_'+stem+'_plot.npz'),**data)
    def map_plot(values,label,limits,stem,category,bounds=None,overlay=None,quiver=None):
        signed=limits[0]<0;levels=np.linspace(*limits,41)
        if signed:levels[20]=0
        fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained')
        m=ax.contourf(x,y,np.ma.masked_invalid(values),levels=levels,norm=Normalize(*limits),cmap='RdBu_r' if signed else 'viridis',corner_mask=False,antialiased=False)
        cb=fig.colorbar(m,ax=ax,label=label);cb.locator=MaxNLocator(nbins=6);cb.update_ticks()
        if stem=='correlation':cb.set_ticks([-1,-.5,0,.5,1])
        ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$y_p\ (\mathrm{mm})$',aspect='equal')
        ax.set_xlim(*(bounds or full_bounds)[:2]);ax.set_ylim(*(bounds or full_bounds)[2:]);ax.margins(0);ax.set_autoscale_on(False)
        if category=='report':
            np.testing.assert_allclose(ax.get_xlim(),[-w,w]);np.testing.assert_allclose(ax.get_ylim(),[low,high])
        if overlay:overlay(ax)
        save(fig,stem,category,'Native-grid rigidly rotated field. '+('Common supported rectangular report crop.' if category=='report' else 'Complete rotated measurement support; outside-domain corners remain unmeasured.'),dict(x_p_mm=x,y_p_mm=y,values=values,color_limits=limits,levels=levels,display_bounds=np.array(bounds or full_bounds),actual_xlim_mm=np.array(ax.get_xlim()),actual_ylim_mm=np.array(ax.get_ylim())))
    for field in labels:
        values=a[field];signed=field in ('mean_v_p','mean_u_p','covariance_u_p_v_p','correlation_u_p_v_p')
        lim=1. if field=='correlation_u_p_v_p' else float(np.nanmax(abs(values)))
        if lim==0:lim=1.
        limits=(-lim,lim) if signed else (0,lim)
        if field=='N':
            map_plot(values,labels[field],limits,stems[field],'diagnostic');continue
        map_plot(values,labels[field],limits,stems[field],'full_support')
        map_plot(values,labels[field],limits,stems[field],'report',[-w,w,low,high])
    # Full original field fit is diagnostic, not the main report map.
    fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained');lim=max(float(np.nanmax(abs(z['mean_v']))),float(np.nanmax(abs(a['mean_v_p']))))
    m=ax.contourf(z['x_mm'],z['y_mm'],np.ma.masked_invalid(z['mean_v']),levels=np.linspace(-lim,lim,41),cmap='RdBu_r',corner_mask=False)
    for selected,col,label in [(False,'0.5','Centres outside fit'),(True,'black','Centres used in fit')]:
        group=[r for r in centres if r['in_primary_fit']==selected];ax.plot([r['x_peak_mm'] for r in group],[r['y_lab_mm'] for r in group],'.',ms=2,color=col,label=label)
    yy=np.array([ymin,ymax]);ax.plot(fit['intercept_x_at_lab_y0_mm']+fit['slope_dx_dy']*yy,yy,'k--',lw=.8,label='One fitted straight axis');ax.plot(*origin,'ko',ms=3)
    ax.set(xlabel=r'$x\ (\mathrm{mm})$',ylabel=r'$y\ (\mathrm{mm})$',aspect='equal');ax.legend(frameon=False,fontsize=7);fig.colorbar(m,ax=ax,label=r'$\overline{v}\ (\mathrm{m\,s^{-1}})$');save(fig,'original_field_axis_fit','diagnostic','Extracted local peak positions and one straight fitted axis; origin marks measurement reference plane, not surface.')
    fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,xx,yy,vv,xlabel,ylabel in [(axes[0],z['x_mm'],z['y_mm'],z['mean_v'],'x','y'),(axes[1],x,y,a['mean_v_p'],'x_p','y_p')]:
        m=ax.contourf(xx,yy,np.ma.masked_invalid(vv),levels=np.linspace(-lim,lim,41),cmap='RdBu_r',corner_mask=False);ax.set(xlabel=rf'${xlabel}\ (\mathrm{{mm}})$',ylabel=rf'${ylabel}\ (\mathrm{{mm}})$',aspect='equal')
    fig.colorbar(m,ax=list(axes),label=r'Axial mean velocity $(\mathrm{m\,s^{-1}})$');save(fig,'original_vs_plume_aligned','diagnostic','Laboratory vertical and plume-aligned axial components are distinct projections of the same vectors.')
    fig,ax=plt.subplots(figsize=(6.4,4.3),layout='constrained');ax.plot(axis['y_p_mm'],axis['mean_v_p_m_s'],color='#d62728',label=r'Fixed fitted axis: $x_p=0$');ax.set(xlabel=r'$y_p\ (\mathrm{mm})$',ylabel=labels['mean_v_p']);ax.legend(frameon=False)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(),axis['mean_v_p_m_s']);save(fig,'fitted_axis_axial_velocity_profile','report','Mean axial velocity along the fixed fitted axis, sampled bilinearly with four finite corners and no extrapolation. This is not the original fixed laboratory-column line.',axis)
    for component,field,stem in [('mean_v_p_m_s','mean_v_p','transverse_axial_velocity_profiles'),('mean_u_p_m_s','mean_u_p','transverse_transverse_velocity_profiles')]:
        fig,ax=plt.subplots(figsize=(6.4,4.4),layout='constrained')
        for section in sections:
            ax.plot(section['x_p_mm'],section[component],label=rf'$y_p={section["y_p_mm"][0]:g}$ mm')
            np.testing.assert_array_equal(ax.lines[-1].get_ydata(),section[component])
        ax.axvline(0,color='0.5',ls=':',lw=.6);ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=labels[field],xlim=(-w,w));ax.legend(frameon=False,fontsize=8,ncol=2);save(fig,stem,'report','Fixed plume-aligned axial stations measured from the measurement reference plane; full-supported section values are exported, including signed components and gaps.')
    fig,ax=plt.subplots(figsize=(6.4,4.3),layout='constrained');ax.plot(fixed_y,original['mean_v'],':',color='.5',label='Fixed laboratory column, laboratory v');ax.plot(fixed_y,fixed_v,'k--',label='Same column, aligned component');ax.plot(axis['y_p_mm'],axis['mean_v_p_m_s'],color='#d62728',label='Fixed fitted axis, aligned component');ax.set(xlabel=r'$y_p\ (\mathrm{mm})$',ylabel=r'Axial velocity $(\mathrm{m\,s^{-1}})$');ax.legend(frameon=False,fontsize=8);save(fig,'different_profile_lines_comparison','diagnostic','Different spatial extraction lines; original laboratory-column exports preserved byte-for-byte under reference/.')
    fig,ax=plt.subplots(figsize=(6.4,4.3),layout='constrained');ax.plot([r['y_lab_mm'] for r in centres],[r['residual_x_mm'] for r in centres],'.-',ms=2,lw=.7);ax.axhline(0,color='k',lw=.6);ax.axvspan(fit['actual_lower_y_mm'],fit['actual_upper_y_mm'],alpha=.1);ax.set(xlabel=r'$y\ (\mathrm{mm})$',ylabel=r'$x_c-x_{\mathrm{fit}}\ (\mathrm{mm})$');save(fig,'straight_axis_residuals','diagnostic','Residual structure is retained; interval sensitivity is not a confidence interval.')
    # Already calibrated actual instantaneous frames, never reconstructed from means.
    frames=[];instant_checks=[]
    for path in sorted((out.parent/'data').glob('frame_*.npz')):
        with np.load(path) as f:frame={k:f[k].copy() for k in f.files}
        np.testing.assert_array_equal(frame['x_mm'],z['x_mm']);np.testing.assert_array_equal(frame['y_mm'],z['y_mm'])
        vel=np.stack([frame['u'],frame['v']],axis=-1);rot=vel@Q.T
        np.testing.assert_allclose(rot@Q,vel,rtol=1e-12,atol=1e-14);np.testing.assert_allclose(np.linalg.norm(rot,axis=-1),np.linalg.norm(vel,axis=-1),rtol=1e-13,atol=1e-14)
        np.testing.assert_array_equal(np.isfinite(rot).all(axis=-1),np.isfinite(vel).all(axis=-1))
        data=dict(x_p_mm=x,y_p_mm=y,u_p=rot[...,0],v_p=rot[...,1],original_lab_x_mm=frame['x_mm'],original_lab_y_mm=frame['y_mm'],original_lab_u=frame['u'],original_lab_v=frame['v'],native_row=np.indices(x.shape)[0],native_column=np.indices(x.shape)[1])
        np.savez_compressed(out/'data'/(path.stem+'_plume_aligned.npz'),**data);frames.append((path.stem,rot));instant_checks.append(dict(source=str(path),sha256=sha(path),inverse_and_magnitude_verified=True))
    shared=[max(float(np.nanmax(abs(v[...,k]))) for _,v in frames) for k in (0,1)]
    ref=float(np.ceil(max(shared)))
    for stem,rot in frames:
        for k,name,symbol in [(0,'transverse','u'),(1,'axial','v')]:
            for category,bounds in [('report',[-w,w,low,high]),('full_support',full_bounds)]:map_plot(rot[...,k],rf'${symbol}_p\ (\mathrm{{m\,s^{{-1}}}})$',(-shared[k],shared[k]),stem+'_'+name+'_velocity',category,bounds)
        for category,bounds in [('report',[-w,w,low,high]),('full_support',full_bounds)]:
            fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained');sl=(slice(None,None,6),slice(None,None,3));valid=np.isfinite(rot).all(axis=-1)
            q=ax.quiver(x[sl],y[sl],np.ma.array(rot[...,0],mask=~valid)[sl],np.ma.array(rot[...,1],mask=~valid)[sl],angles='xy',scale_units='xy',scale=ref/(.07*max(np.ptp(x),np.ptp(y))),width=.003)
            ax.quiverkey(q,.75,.94,ref,rf'${ref:g}\ \mathrm{{m\,s^{{-1}}}}$',coordinates='axes',labelpos='S');ax.set(xlabel=r'$x_p\ (\mathrm{mm})$',ylabel=r'$y_p\ (\mathrm{mm})$',aspect='equal',xlim=bounds[:2],ylim=bounds[2:]);save(fig,stem+'_quiver',category,'Actual instantaneous rotated components; shared physical arrow scale across representative frames. Plot subsampling changes arrows displayed only, not exported native data.')
    # Complete numerical table includes tensor quantities and inverse-mapped source coordinates.
    rows,cols=np.indices(x.shape)
    fields=dict(native_row=rows.ravel(),native_column=cols.ravel())
    for key,value in a.items():
        if np.shape(value)==x.shape:fields[key]=value.ravel()
    csv_out(out/'data/plume_aligned_mean_fields.csv',fields)
    json_out(out/'numerical_manifest.json',dict(files=[dict(path=str(p.relative_to(out)).replace('\\','/'),sha256=sha(p)) for p in sorted((out/'data').glob('*')) if p.is_file()]))
    json_out(out/'validation/complete_render_checks.json',dict(report_crop_corners_inside_domain=True,profile_plot_arrays_verified=True,instantaneous_frames=instant_checks,shared_instantaneous_limits=[[-v,v] for v in shared],no_extrapolation=True,internal_masks_retained=True))
    return manifest


def write_complete_report(out,cfg,m,stations,manifest,checks):
    from .plume_alignment import json_out
    fit=m['primary_fit'];crop=m['report_crop'];origin=m['reference_point_lab_mm'];angles=[s['angle_right_of_lab_vertical_deg'] for s in m['fit_sensitivity']]
    command=r'.\.venv\Scripts\python.exe -m PZT_Atomiser_Analysis.analysis.update_plume_aligned_notebook'
    text=f'''# Complete plume-aligned measurement package

Interaction predictions remain paused. No crossing location is selected or updated.

## Reproduction
From the local PZT_droplets project root:
```powershell
{command}
```
This executes the alignment notebook section, rebuilds all outputs and saves the executed notebooks. Existing generated outputs are archived under Plume_Aligned/Archive/<timestamp> before replacement. Configuration: config/plume_alignment.json. Source: analysis/plume_alignment.py and analysis/plume_alignment_complete.py.

## Measurements and axis estimate
Already calibrated authoritative laboratory fields are used unchanged: 5.35 micrometres/pixel, 25 microseconds, 0.214 m/s per pixel/pair. No calibration is applied here. Original finite-pair sample counts and masks are preserved. Original images linked to this PIV sequence were unavailable; previously inspected nearby EVOS microscope images and assembly illustration do not establish acquisition tilt. The archived image_inventory.json documents that inspection.
One subgrid quadratic peak of laboratory mean_v is extracted per resolved native y row; equal-weight least squares fits x_c=a+b*y. Fit angle: {fit['angle_right_of_lab_vertical_deg']:.9f} degrees right of laboratory vertical. Actual fit interval: {fit['actual_lower_y_mm']:.4f} to {fit['actual_upper_y_mm']:.4f} mm, {fit['centre_count']} centres. RMS horizontal residual {fit['residual_RMS_x_mm']:.7f} mm; maximum absolute residual {fit['residual_max_abs_x_mm']:.7f} mm. Interval sensitivity: {min(angles):.6f} to {max(angles):.6f} degrees. This is a hydrodynamic straight-axis estimate, not independently established mounting tilt. Residual curvature and asymmetry are not attributed entirely to mounting. Subgrid localisation does not imply corresponding instrument accuracy. Secondary half-maximum centroid fitting and interval results are recorded in transformation.json.

## One rigid transformation and measurement reference
Origin in laboratory coordinates: ({origin[0]:.9f}, {origin[1]:.9f}) mm, the fitted-axis intersection with the lowest supported laboratory-y measurement plane. It is the measurement reference plane, not the physical atomiser surface. x_p is transverse, y_p increases along the plume; y_p is not laboratory vertical.
Position p_p=Q(p_lab-p_ref); velocity w_p=Q*w_lab; covariance C_p=Q*C_lab*Q^T. Inverse p_lab=Q^T*p_p+p_ref and w_lab=Q^T*w_p. Q and its inverse, original lab coordinates, calibrated lab source paths and transformed gravity {m['gravity_plume_m_s2']} m/s^2 are exported. RMS values are square roots of rotated tensor diagonal entries, never rotated as a vector. No rows are shifted independently, mirrored, averaged or symmetrised; nonzero transverse velocity and actual curvature remain. This representation does not establish what physical remounting would produce.

## Native support and report window
Principal standalone corrected maps, all using identical actual rectangular plot limits: figures/report/. Uncropped native-grid plots are diagnostic references only: diagnostics/full_support/. Before/after, fit, residuals, sample counts and original-column comparisons: diagnostics/figures/.
Report x_p bounds {crop['x_p_bounds_mm']} mm; y_p bounds {crop['y_p_bounds_mm']} mm. Width covers all extracted half-maximum lobes plus the configured margin, rounded upward; given this width, the longest axial interval with all four corners in measured laboratory support is retained. Retained axial span {100*crop['retained_fitted_axis_fraction']:.3f}% of the fitted-axis measured span and {100*crop['retained_axial_fraction']:.3f}% of full rotated bounding-box vertical extent; excluded {crop['excluded_axial_mm']:.6f} mm. Native nodes inside the report rectangle: {100*crop['native_node_fraction_inside_crop']:.2f}%. This node fraction includes deliberate lateral crop, not data rejection. All native nodes remain exported and shown in full-support maps. Report windows clip contour polygons; no window stretching, missing-corner filling or extrapolation occurs. Internal missing-data masks are preserved (corner_mask=False). Spatial aspect ratio is equal and transverse display bounds symmetric. Signed maps use symmetric linear blue-white-red normalization with uniformly spaced levels; matching full/report maps share limits, instantaneous components share limits across all selected frames. RMS starts at zero; correlation spans [-1,1].

## Profiles and instantaneous measurements
The fitted-axis profile samples the fixed straight x_p=0 line at native laboratory y planes, using four-finite-corner bilinear interpolation without extrapolation. It differs from the original fixed laboratory-column line, whose original exports are copied byte-for-byte to reference/. Transverse cuts are fixed y_p stations {cfg['transverse_stations_plume_y_mm']} mm; source lab positions, signed components and finite support are saved. Full section values are exported even where the report display clips lateral tails. No smoothing or adjacent-column averaging is used.
Selected saved instantaneous frames are discovered from original data/frame_*.npz (currently 1, 879, 1757). Actual already-calibrated u/v arrays are rotated, not reconstructed from means. Both component maps and separate quivers are exported in full support and report form. Quivers subsample displayed arrows only and share a physical scale.

## Verification and retained asymmetry
Inverse position/vector transforms, vector-magnitude preservation, covariance inverse/trace/eigenvalues, unchanged native N and masks, outside-support sampling, and direct profile-plot/export agreement are checked. Original laboratory numerical data and figures are protected by SHA256 comparisons. Instantaneous inverse/magnitude/mask checks are recorded separately. No interaction integration is run.
Station metrics in data/transverse_profile_metrics.csv retain left/right half-width differences; genuine asymmetry and nonzero transverse velocity remain. These descriptive differences are not a formal statistical significance test.

## Outputs
{len(manifest)} separate figures, each PDF, SVG and 600-dpi PNG. report_manifest.json enumerates exact paths and captions; numerical native grids, original lab coordinates, profiles, full tensors, instantaneous arrays, support and transformations are in data/ and transformation.json.
'''
    (out/'README.md').write_text(text,encoding='utf-8')
    latex=r'''\section{Plume-aligned measured PZT field}
The calibrated laboratory measurements are represented by a single rigid transformation,
$\mathbf{x}_p=Q(\mathbf{x}-\mathbf{x}_{\rm ref})$ and $\mathbf{v}_p=Q\mathbf{v}$,
with $C_p=QCQ^T$. The measurement reference is the fitted-axis intersection with
the lowest supported laboratory height, not the atomiser surface. The axial coordinate
$y_p$ is plume-aligned, not laboratory vertical. Native coordinates and inverse mapping
are retained; contours are not resampled. Line samples use finite-corner bilinear
interpolation without extrapolation. The fitted-axis line is distinct from the original
fixed laboratory-column extraction. Curvature, asymmetry and nonzero transverse velocity
remain. This coordinate transformation does not predict physical remounting.
'''
    latex+=f'Fitted angle: {fit["angle_right_of_lab_vertical_deg"]:.6f} degrees; fit interval {fit["actual_lower_y_mm"]:.4f}--{fit["actual_upper_y_mm"]:.4f} mm; horizontal residual RMS {fit["residual_RMS_x_mm"]:.6f} mm.\n'
    latex+=f'Laboratory reference: ({origin[0]:.9f}, {origin[1]:.9f}) mm. Report crop: $x_p\\in[{crop["x_p_bounds_mm"][0]:g},{crop["x_p_bounds_mm"][1]:g}]$ mm and $y_p\\in[{crop["y_p_bounds_mm"][0]:.6f},{crop["y_p_bounds_mm"][1]:.6f}]$ mm. The common crop maximises axial support for the selected symmetric width; full native support is exported separately.\n'
    latex+='\nAlready calibrated data use 5.35 micrometres/pixel and 25 microseconds; calibration is not reapplied. Interaction predictions remain paused.\n'
    (out/'methods.tex').write_text(latex,encoding='utf-8')
    lines=[]
    for item in manifest:
        quantity=item['stem'].replace('_',' ')
        caption=quantity[0].upper()+quantity[1:]+'. '+item['caption']+' Transverse $x_p$ and plume-aligned axial $y_p$ are referenced to the lowest measured laboratory-y plane, not the atomiser surface. One rigid transformation preserves genuine asymmetry.'
        lines.append('% '+item['path']+'\n\\caption{'+caption+'}\n')
    (out/'captions.tex').write_text('\n'.join(lines),encoding='utf-8')
