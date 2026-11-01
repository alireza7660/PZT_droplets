"""Validate actual proposed cases and geometry-independent dimensionless outputs."""
import json
import numpy as np
from PIL import Image
from .combined_flow_prediction import ROOT,PROJECT,sources_and_arrays,write_json,sha
from .droplet_dynamics import relaxation_time,reynolds
from .prescribed_gas import MeasuredGasField


def read(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))


def main():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    pzt,jet,xj,yj,D0,source=sources_and_arrays(cfg);props=json.loads((out/'fluid_properties.json').read_text())
    d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,props))
    dim=out/'Dimensionless_Numbers';study=out/'Design_Study'
    z=read(dim/'znmf_local_centreline_numbers.csv');good=z['definition_supported'].astype(bool)
    parameters=json.loads((dim/'parameters_and_validation.json').read_text());f=parameters['actuation_frequency_Hz'];nu=parameters['nu_air_m2_s']
    np.testing.assert_allclose(z['Re_D_c'][good],z['Uc_measured_m_s'][good]*D0/nu,rtol=1e-14)
    np.testing.assert_allclose(z['Sr_D_c'][good],f*D0/z['Uc_measured_m_s'][good],rtol=1e-14)
    np.testing.assert_allclose(z['Stk_D'][good],tau*z['Uc_measured_m_s'][good]/D0,rtol=1e-14)
    np.testing.assert_allclose(z['Stk_f'],tau*f,rtol=1e-14)
    assert np.isnan(z['Re_D_c'][~good]).all() and np.isnan(z['Sr_D_c'][~good]).all()
    rp=read(dim/'pzt_zero_background_reference_field.csv')
    np.testing.assert_array_equal(rp['mean_u_m_s'],pzt['mean_u'].ravel());np.testing.assert_array_equal(rp['mean_v_m_s'],pzt['mean_v'].ravel())
    np.testing.assert_allclose(rp['Re_p_zero_gas_reference'],reynolds(d,np.c_[rp['mean_u_m_s'],rp['mean_v_m_s']],[0,0],props),rtol=1e-14)
    cases=read(study/'case_comparison.csv');placements=read(study/'placement_comparison.csv')
    count=0;terminations={}
    for case in cases:
        q=read(study/(case['case']+'.csv'));i,j=int(case['release_row']),int(case['release_column'])
        np.testing.assert_allclose([q['droplet_u_m_s'][0],q['droplet_v_m_s'][0]],[pzt['mean_u'][i,j],pzt['mean_v'][i,j]],rtol=0,atol=1e-15)
        np.testing.assert_allclose([q['x_mm'][0],q['y_mm'][0]],[pzt['x_mm'][i,j],pzt['y_mm'][i,j]],rtol=0,atol=1e-12)
        assert q['t_s'][0]==0 and q['t_s'][-1]>0 and len(q)>1
        field=MeasuredGasField(xj,yj,jet['U'],jet['V'],np.array([case['nozzle_x_mm'],case['nozzle_y_mm']])*1e-3,0)
        gas=np.array([field.sample(pos*1e-3) for pos in np.c_[q['x_mm'],q['y_mm']]])
        np.testing.assert_allclose(gas,np.c_[q['gas_u_m_s'],q['gas_v_m_s']],rtol=1e-12,atol=1e-12)
        drop=np.c_[q['droplet_u_m_s'],q['droplet_v_m_s']]
        np.testing.assert_allclose(q['Re_in_plane'],reynolds(d,drop,gas,props),rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(q['St_D_local'],tau*np.linalg.norm(gas,axis=1)/D0,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(q['horizontal_deflection_mm'],q['x_mm']-q['x_mm'][0],atol=1e-12)
        np.testing.assert_allclose(q['vertical_displacement_mm'],q['y_mm']-q['y_mm'][0],atol=1e-12)
        np.testing.assert_allclose(q['tau_stokes_s'],tau,rtol=1e-14)
        terminations[case['termination']]=terminations.get(case['termination'],0)+1;count+=1
    assert count==len(placements)*len(cfg['design_study']['release_column_offsets'])
    for check in json.loads((study/'convergence.json').read_text()):
        errors=check['max_errors_x_y_u_v'];assert max(errors[:2])<cfg['solver']['position_convergence_tolerance_m'];assert max(errors[2:])<cfg['solver']['velocity_convergence_tolerance_m_s']
    file_lists={}
    for folder in (dim,study):
        manifest=json.loads((folder/'manifest.json').read_text())
        for stem in manifest['figures']:
            for ext in ('pdf','png'):assert (folder/(stem+'.'+ext)).stat().st_size>1000
            with Image.open(folder/(stem+'.png')) as im:assert min(im.info['dpi'])>599
        file_lists[folder.name]=sorted(p.name for p in folder.iterdir() if p.is_file())
    result=dict(successfully_integrated_cases=count,proposed_placements=len(placements),termination_counts=terminations,
        release_measured_components_and_positions_verified=True,all_trajectory_gas_points_supported=True,
        evolving_slip_Re_and_local_Stokes_verified=True,half_step_convergence_passed=True,
        dimensionless_positive_measured_gaps_preserved=True,actuation_frequency_Hz=f,Stk_f=tau*f,
        zero_background_reference_not_measured_gas=True,original_PZT_signed_components_preserved=True)
    write_json(out/'validation/prospective_studies_verification.json',result)
    write_json(out/'prospective_generated_files.json',file_lists)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
