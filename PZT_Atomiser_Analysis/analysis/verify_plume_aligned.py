"""Independent saved-output checks for the complete plume-aligned package."""
import json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
from .plume_alignment import ROOT,PROJECT,sha,json_out,sample_native


def verify_outputs(out=None):
    out=Path(out or PROJECT/'Final_Results/Batch4_Cam1/Plume_Aligned')
    metadata=json.loads((out/'transformation.json').read_text());figures=json.loads((out/'report_manifest.json').read_text())['figures']
    source=np.load(out.parent/'data/statistics.npz');a=np.load(out/'data/statistics_plume_aligned.npz');Q=np.array(metadata['orthonormal_matrix_lab_to_plume']);origin=np.array(metadata['reference_point_lab_mm'])
    np.testing.assert_array_equal(a['N'],source['N']);assert a['N'].shape==source['N'].shape
    np.testing.assert_allclose(np.stack([a['x_p_mm'],a['y_p_mm']],-1)@Q+origin,np.stack([source['x_mm'],source['y_mm']],-1),atol=1e-13)
    native_y=source['y_mm'][np.isfinite(source['mean_u'])&np.isfinite(source['mean_v'])].min();assert origin[1]==native_y
    profile=np.load(out/'data/fitted_axis_axial_profile.npz');values=np.genfromtxt(out/'data/fitted_axis_axial_profile.csv',delimiter=',',names=True)
    np.testing.assert_allclose(profile['x_p_mm'],0,atol=1e-12)
    for name in ('y_p_mm','lab_x_mm','lab_y_mm','mean_u_p_m_s','mean_v_p_m_s'):np.testing.assert_allclose(profile[name],values[name],atol=1e-14)
    sampled=sample_native(source,np.c_[profile['lab_x_mm'],profile['lab_y_mm']])@Q.T
    np.testing.assert_allclose(sampled[:,1],profile['mean_v_p_m_s'],atol=1e-14)
    plot=np.load(out/'data/report_fitted_axis_axial_velocity_profile_plot.npz')
    np.testing.assert_array_equal(plot['mean_v_p_m_s'],profile['mean_v_p_m_s'])
    # Explicit internal-gap and outside-support tests; these operate on disposable copies.
    copy={k:source[k].copy() for k in source.files};ii,jj=np.array(copy['mean_u'].shape)//2;copy['mean_u'][ii,jj]=np.nan
    point=[[copy['x_mm'][ii,jj],copy['y_mm'][ii,jj]]];assert np.isnan(sample_native(copy,point)).all()
    assert np.isnan(sample_native(source,[[source['x_mm'].min()-1,source['y_mm'].min()-1]])).all()
    bounds=np.array(metadata['report_crop']['x_p_bounds_mm']+metadata['report_crop']['y_p_bounds_mm'])
    checked=0
    for path in (out/'data').glob('report_*_plot.npz'):
        saved=np.load(path)
        if 'actual_xlim_mm' not in saved.files:continue
        np.testing.assert_allclose(saved['actual_xlim_mm'],bounds[:2]);np.testing.assert_allclose(saved['actual_ylim_mm'],bounds[2:]);checked+=1
    assert checked>=6
    # The complete rectangle perimeter is within measured lab support, not a rotated bounding box.
    xx=np.linspace(bounds[0],bounds[1],100);yy=np.linspace(bounds[2],bounds[3],200)
    edge=np.vstack([np.c_[xx,np.full_like(xx,bounds[2])],np.c_[xx,np.full_like(xx,bounds[3])],np.c_[np.full_like(yy,bounds[0]),yy],np.c_[np.full_like(yy,bounds[1]),yy]])
    lab=edge@Q+origin
    assert lab[:,0].min()>=source['x_mm'].min()-1e-10 and lab[:,0].max()<=source['x_mm'].max()+1e-10
    assert lab[:,1].min()>=source['y_mm'].min()-1e-10 and lab[:,1].max()<=source['y_mm'].max()+1e-10
    for stem in ('mean_axial_velocity','mean_transverse_velocity','axial_RMS','transverse_RMS','axial_transverse_covariance','correlation'):
        full=np.load(out/'data'/('full_support_'+stem+'_plot.npz'));report=np.load(out/'data'/('report_'+stem+'_plot.npz'))
        np.testing.assert_array_equal(full['values'],report['values']);np.testing.assert_array_equal(full['color_limits'],report['color_limits'])
    for ext in ('csv','npz'):assert sha(out/'reference'/('original_fixed_column_profile.'+ext))==sha(out.parent/'data'/('mean_v_vertical_profile.'+ext))
    protected=json.loads((out/'validation/protected_file_hashes.json').read_text());assert all(sha(p)==h for p,h in protected.items())
    for item in figures:
        p=out/item['path']
        for ext in ('.pdf','.svg','.png'):assert p.with_suffix(ext).stat().st_size>1000
        with Image.open(p.with_suffix('.png')) as image:assert min(image.info['dpi'])>599
    paths=[(out/i['path']).with_suffix('.png') for i in figures]
    for batch in range(0,len(paths),12):
        canvas=Image.new('RGB',(1600,1350),'white');draw=ImageDraw.Draw(canvas)
        for i,p in enumerate(paths[batch:batch+12]):
            with Image.open(p) as image:
                image.thumbnail((390,420));xx=(i%4)*400;yy=(i//4)*450;canvas.paste(image,(xx,yy+24));draw.text((xx+3,yy+3),p.parent.name+'/'+p.stem,fill='black')
        canvas.save(out/'validation'/f'visual_review_{batch//12+1}.png')
    result=dict(rectangle_actual_limits_verified_in_contours=checked,rectangle_entire_perimeter_in_measured_domain=True,figure_count=len(figures),report_figures=sum(i['category']=='report' for i in figures),full_support_figures=sum(i['category']=='full_support' for i in figures),diagnostic_figures=sum(i['category']=='diagnostic' for i in figures),all_PDF_SVG_600dpi_PNG_exports_verified=True,protected_files_unchanged=len(protected),measurement_plane_origin_verified=True,axis_profile_CSV_NPZ_source_and_plot_agree=True,finite_profile_sample_count=int(np.isfinite(profile['mean_v_p_m_s']).sum()),internal_gap_and_outside_support_tests_pass=True,report_full_support_limits_identical=True,original_column_profile_byte_identical=True,interaction_predictions_run=False)
    json_out(out/'validation/saved_output_verification.json',result)
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':verify_outputs()
