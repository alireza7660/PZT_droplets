"""Active co-directed prediction from full-support plume-aligned measurements."""
from .aligned_codirected import run,load_context,prepare


def notebook_cells():
    import nbformat as nb
    first=nb.v4.new_markdown_cell('## Co-directed response using corrected plume-aligned measurements\n\nThe original laboratory-height/source-distance pairing identifies a supported overtaking crossing and a subsequent native entry plane. Both aligned droplet components initialise already-developed droplets there. One rigid jet placement, transformed gravity, finite-Re drag and strict measured support determine the subsequent response. The diameter is a nominal 10 micrometre model estimate, not measured sizing.')
    code=nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.co_directed import run as run_codirected
prediction=run_codirected()
print(json.dumps(prediction,indent=2))
package=output/'Combined_Flow_Prediction/Co_Directed'
for name in ['source_comparison_reproduced','measured_and_predicted_axial_comparison','trajectories_upward_jet','trajectories_entry_detail']:
    display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    text=nb.v4.new_markdown_cell('## Slip, response and local jet dimensionless numbers\n\nAll trajectory comparisons use the same supported ascending axial interval. Re_p uses vector slip; Stk_D,local uses tau_0 and D0 divided by the local gas-vector magnitude. Sr denotes Strouhal number, not Stokes number. Jet forcing is 120 Hz, distinct from 10-Hz PIV acquisition. Measured PZT-only means are an ensemble baseline, not individual droplet histories.')
    plots=nb.v4.new_code_cell('''print((package/'data/initial_conditions.csv').read_text())
print((package/'data/common_height_comparison.csv').read_text())
print((package/'data/scalar_quantities.csv').read_text())
for name in ['axial_slip','particle_Re','particle_Stokes','axial_acceleration','jet_local_Reynolds','jet_local_Strouhal']:
    display(Image(filename=str(package/'figures'/(name+'.png')),width=700))
''')
    for cell in (first,code,text,plots):cell.metadata['co_directed']=True;cell.metadata['combined_flow_prediction']=True
    return [first,code,text,plots]
