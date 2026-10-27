"""Create and execute the Cartesian, corrected-calibration final notebook."""
import asyncio
import json
import os
from pathlib import Path
import shutil
import sys
import nbformat as nb
from nbclient import NotebookClient
from jupyter_client import AsyncKernelManager
from jupyter_client.kernelspec import KernelSpecManager


def main():
    root=Path(__file__).resolve().parents[2]
    project=root/'PZT_Atomiser_Analysis'
    runtime=project/'Final_Results/Batch4_Cam1/diagnostics/jupyter_runtime'
    runtime.mkdir(parents=True,exist_ok=True);os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    b=nb.v4.new_notebook()
    active=json.loads((project/'config/analysis.json').read_text())
    b.metadata.kernelspec={'display_name':'Python (PZT workspace)','language':'python','name':'python3'}
    b.cells=[nb.v4.new_markdown_cell(f'''# PZT atomiser final analysis

The PZT calibration is {active['confirmed_spatial_calibration_um_per_pixel']:g} micrometres/pixel with pulse separation
{active['confirmed_pulse_separation_seconds']*1e6:g} microseconds: instantaneous pixel displacements are converted once by
{active['velocity_scale']:.3f} m/s. This calibration supersedes both 5.05 and 7.25.
Cartesian x/u is horizontal/transverse; y/v is vertical/axial.
Naming changes do not transpose, rotate or reflect the field. The exported
image-reference origin and existing station indices/frame keys are retained.

Mean and instantaneous axial contours, including the station map, display
all finite signed measurements. Blue is negative, white is zero and red is
positive. Symmetric limits and one linear normalization are used, with
uniformly spaced contour levels and no separate stretching of signs. Each
mean component uses its own maximum absolute ensemble value, excluding flag
diagnostics. Selected instantaneous fields share limits. Only missing or
invalid data are masked. Values are recalculated from the original displacement
samples; sample selection, station indices, averaging and orientation are retained.
'''),nb.v4.new_code_cell('''from pathlib import Path
import sys,json
root=Path.cwd().resolve()
while not (root/'RAW_data').exists():
    if root==root.parent: raise RuntimeError('Workspace root not found')
    root=root.parent
sys.path.insert(0,str(root))
config_path=root/'PZT_Atomiser_Analysis/config/analysis.json'
config=json.loads(config_path.read_text())
from math import isclose
assert isclose(config['velocity_scale'],config['confirmed_spatial_calibration_um_per_pixel']*1e-6/config['confirmed_pulse_separation_seconds'])
print(json.dumps(config,indent=2))
'''),nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.pipeline import run
results=run(config_path)
output=root/config['output']
'''),nb.v4.new_code_cell('''from IPython.display import display,Image,Markdown
display(Markdown((output/'methods/METHODS.md').read_text(encoding='utf-8')))
display(Markdown('## Primary mean-component contours'))
for name in ['mean_v','mean_u']:
    display(Image(filename=str(output/'figures'/(name+'.png')),width=600))
display(Markdown('## Component fluctuations, profiles and distributions'))
for name in ['rms_v','rms_u','covariance','correlation','sampling_stations',
             'profile_mean_v','profile_rms_v','profile_rms_u','profile_covariance',
             'v_pdf','u_pdf','v_cdf','u_cdf','frame_879_v','frame_879_u','frame_879_quiver']:
    print(name)
    display(Image(filename=str(output/'figures'/(name+'.png')),width=600))
display(Markdown('## Supplementary in-plane velocity magnitude'))
for name in ['magnitude_mean_vector','mean_instantaneous_velocity_magnitude','velocity_magnitude_pdf']:
    print(name)
    display(Image(filename=str(output/'figures'/(name+'.png')),width=600))
'''),nb.v4.new_code_cell('''print('Measured extrema and display limits are reported separately:')
print((output/'methods/summary.json').read_text())
print((output/'methods/validation.json').read_text())
import csv
with (output/'data/vector_comparison.csv').open() as f: rows=list(csv.DictReader(f))
columns=list(rows[0])
table=[' | '.join(columns),' | '.join(['---']*len(columns))]
table+=[' | '.join(row[c] for c in columns) for row in rows]
display(Markdown(chr(10).join(table)))
'''),nb.v4.new_markdown_cell('''## Reuse

Edit config/analysis.json and rerun. The same command-line entry point calls
the same reusable accumulator and plotting functions. The eight signed
statistics and numerical figure sources are saved separately from their
display masks. See methods/figure_manifest.json for every export/data link,
methods/colour_limits.json for shared limits, and methods/CAPTIONS.md or .tex
for captions. Previous calibrated archives and original HDF5 files are retained.
''')]
    target=project/'notebooks/PZT_final_analysis.ipynb';target.parent.mkdir(exist_ok=True)
    nb.write(b,target)
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'];manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(b,timeout=1200,resources={'metadata':{'path':str(root)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager)
    try:client.execute()
    finally:nb.write(b,target)
    shutil.copy2(target,project/'notebooks/Spatial_maps.ipynb')
    shutil.copy2(target,project/'Final_Results/Batch4_Cam1/methods'/target.name)
    print('Executed final notebook:',target)


if __name__=='__main__':main()
