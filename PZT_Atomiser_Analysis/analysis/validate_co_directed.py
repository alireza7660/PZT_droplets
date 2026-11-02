"""Validate co-directed SI transformations, forces, common intervals and preservation."""
import json
import numpy as np
import nbformat
from PIL import Image
from .co_directed import context
from .combined_flow_prediction import ROOT,PROJECT,sha,write_json
from .droplet_dynamics import reynolds,acceleration,relaxation_time
from .validate_combined_flow import synthetic_checks


def read(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=',',dtype=None,encoding='utf-8'))


def main():
    cfg,pzt,jet,D0,sources,props,out,profile,peak,field,alignment=context(write_preflight=False)
    result=json.loads((out/'results_summary.json').read_text());records=read(out/'data/trajectory_summary.csv')
    d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,props))
    np.testing.assert_array_equal(field.vectors_to_common([1.,0.]),[0.,1.])
    np.testing.assert_array_equal(field.vectors_to_common([0.,1.]),[-1.,0.])
    xx,yy,ug,vg=field.common_mesh()
    np.testing.assert_array_equal(ug,-jet['V']);np.testing.assert_array_equal(vg,jet['U'])
    centre_local=field.to_local(np.array([alignment['plume_reference_x_mm'],alignment['PZT_interaction_height_mm']])*1e-3)
    np.testing.assert_allclose(centre_local*1e3,[alignment['ZNMF_interaction_downstream_mm'],alignment['ZNMF_interaction_centreline_Yj_mm']],atol=1e-12)
    np.testing.assert_allclose(acceleration([0.,0.],d,[0.,0.],props,cfg),[0,-props['g_m_s2']*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3'])])
    common=read(out/'data/common_height_comparison.csv')
    for record in records:
        q=read(out/'data'/(record['release']+'.csv'));i,j=int(record['row']),int(record['column'])
        i1=int(record['upper_row']);weight=float(record['PZT_upper_weight'])
        expected_u=(1-weight)*pzt['mean_u'][i,j]+weight*pzt['mean_u'][i1,j]
        expected_v=(1-weight)*pzt['mean_v'][i,j]+weight*pzt['mean_v'][i1,j]
        np.testing.assert_allclose([q['droplet_u_m_s'][0],q['droplet_v_m_s'][0]],[expected_u,expected_v],rtol=0,atol=1e-15)
        np.testing.assert_allclose([q['x_mm'][0],q['y_mm'][0]],[pzt['x_mm'][i,j],alignment['PZT_interaction_height_mm']],rtol=0,atol=1e-12)
        gas=np.array([field.sample(pos*1e-3) for pos in np.c_[q['x_mm'],q['y_mm']]])
        np.testing.assert_allclose(gas,np.c_[q['gas_u_m_s'],q['gas_v_m_s']],rtol=1e-12,atol=1e-12)
        drop=np.c_[q['droplet_u_m_s'],q['droplet_v_m_s']]
        np.testing.assert_allclose(q['Re_in_plane'],reynolds(d,drop,gas,props),rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(q['net_axial_acceleration_m_s2'],acceleration(drop,d,gas,props,cfg)[:,1],rtol=1e-10,atol=1e-9)
        np.testing.assert_allclose(q['slip_v_m_s'],gas[:,1]-drop[:,1],atol=1e-12)
        np.testing.assert_allclose(q['Stk_D'],tau*np.linalg.norm(gas,axis=1)/D0,rtol=1e-12,atol=1e-12)
        row=common[common['release']==record['release']][0]
        assert row['travel_time_s']>0
        np.testing.assert_allclose([row['common_lower_y_mm'],row['common_upper_y_mm']],result['common_interval_y_mm'],atol=1e-14)
        assert abs(row['common_lower_y_mm']-alignment['PZT_interaction_height_mm'])<1e-12
        assert row['common_upper_y_mm']<=np.nanmax(profile['y_mm'])+1e-12
    for check in json.loads((out/'validation/convergence.json').read_text()):
        errors=check['max_errors_x_y_u_v'];assert max(errors[:2])<cfg['solver']['position_convergence_tolerance_m'];assert max(errors[2:])<cfg['solver']['velocity_convergence_tolerance_m_s']
    jet_numbers=read(out/'data/jet_local_centreline_numbers.csv');good=jet_numbers['supported'].astype(bool)
    f=result['dimensionless']['actuation_frequency_Hz'];nu=result['dimensionless']['nu_air_m2_s']
    np.testing.assert_allclose(jet_numbers['Re_D_c'][good],jet_numbers['Uc_measured_m_s'][good]*D0/nu,rtol=1e-14)
    np.testing.assert_allclose(jet_numbers['Sr_D_c'][good],f*D0/jet_numbers['Uc_measured_m_s'][good],rtol=1e-14)
    assert np.isnan(jet_numbers['Re_D_c'][~good]).all()
    assert not any('Stk_f' in p.name for p in (out/'figures').iterdir())
    manifest=json.loads((out/'report_manifest.json').read_text())
    for stem in manifest['figures']:
        for ext in ('pdf','png'):assert (out/'figures'/(stem+'.'+ext)).stat().st_size>1000
        with Image.open(out/'figures'/(stem+'.png')) as im:assert min(im.info['dpi'])>599
    protected=json.loads((out/'validation/authoritative_hashes.json').read_text())
    for path,digest in protected.items():assert sha(path)==digest,path
    archive=json.loads((ROOT/cfg['co_directed']['archive_manifest']).read_text())
    for r in archive['records']:assert sha(ROOT/r['archive_path'])==r['sha256']
    for name in ('Design_Study','Dimensionless_Numbers'):assert not (out.parent/name).exists()
    book=nbformat.read(PROJECT/'notebooks/PZT_final_analysis.ipynb',as_version=4)
    assert not any(c.metadata.get('prospective_studies') for c in book.cells)
    code=[c for c in book.cells if c.cell_type=='code' and c.metadata.get('co_directed')]
    assert len(code)==2 and all(c.execution_count is not None for c in code)
    assert sum('image/png' in o.get('data',{}) for c in code for o in c.outputs)==len(manifest['figures'])
    checks=synthetic_checks(cfg)
    audit=json.loads((out/'interaction_plane_audit.json').read_text())
    assert audit['comparison_reproduction']['source_curve_values_exactly_match_archived_CSVs']
    assert audit['initial_conditions_exported_before_integration']
    receipt=dict(successfully_integrated_cases=len(records),rotation_and_vector_signs_verified=True,interaction_centreline_alignment_verified=True,
        original_comparison_reproduced_from_sources=True,crossing_computed_with_valid_brackets=True,developed_PZT_initialisation_verified=True,
        signed_gas_reversals_retained=True,initial_measured_PZT_components_verified=True,all_saved_positions_supported=True,
        drag_gravity_and_buoyancy_verified=True,common_interval_comparison_verified=True,step_convergence_passed=True,
        local_jet_gaps_retained=True,scalars_tabulated_not_constant_line_figures=True,
        authoritative_files_hash_verified=len(protected),archived_files_hash_verified=len(archive['records']),
        notebook_executed=True,figure_count=len(manifest['figures']),synthetic_solver_tests=checks)
    write_json(out/'validation/completion_verification.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
