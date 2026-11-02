"""Active co-directed prediction: jet target station placed at measured plume base."""
from .aligned_codirected import run,load_context,prepare


def notebook_cells():
    import nbformat as nb
    markdown=nb.v4.new_markdown_cell('## Jet underneath the plume: measured-base entrance\n\nThe first supported fitted-axis measurement defines h=0. The developing ZNMF centreline station Uc=0.96 m/s is placed there using one rigid map Xj=x_j0+h. Both measured near-base droplet components initialise supported starts. The former profile-crossing placement is superseded; positive acceleration is not imposed. The diameter is a nominal 10-micrometre estimate.')
    code=nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.co_directed import run as run_codirected
prediction=run_codirected()
print(json.dumps(prediction,indent=2))
package=output/'Combined_Flow_Prediction/Co_Directed'
first_figures=['source_comparison_reproduced','measured_and_predicted_axial_comparison','reference_entry_response','trajectories_upward_jet','trajectories_entry_detail']
for name in first_figures:
    display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    text=nb.v4.new_markdown_cell('## Force response, minimum, support and dimensionless diagnostics\n\nInitially the reference droplet is faster than the local gas: its axial drag is negative. Zero axial slip and zero net acceleration are separate events; vector slip is recorded. The original seven offsets are audited at exactly h=0, and unsupported starts are rejected without extrapolation. PZT-only means remain a comparison. Stokes and Strouhal numbers have distinct, explicit timescales.')
    figures=nb.v4.new_code_cell('''print((package/'data/entry_mapping.csv').read_text())
print((package/'data/release_candidate_audit.csv').read_text())
print((package/'data/common_height_comparison.csv').read_text())
print((package/'data/reference_event_summary.json').read_text())
for name in json.loads((package/'report_manifest.json').read_text())['figures']:
    if name not in first_figures:
        display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    for cell in (markdown,code,text,figures):cell.metadata['co_directed']=True;cell.metadata['combined_flow_prediction']=True
    return [markdown,code,text,figures]
