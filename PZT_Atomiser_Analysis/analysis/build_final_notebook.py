"""Create and execute the Cartesian, corrected-calibration final notebook."""
import asyncio
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
    b.metadata.kernelspec={'display_name':'Python (PZT workspace)','language':'python','name':'python3'}
    b.cells=[nb.v4.new_markdown_cell('''# PZT atomiser final analysis

The corrected PZT calibration is 7.25 micrometres/pixel with pulse separation
25 microseconds: instantaneous pixel displacements are converted once by
0.290 m/s. Cartesian x/u is horizontal/transverse; y/v is vertical/axial.
Naming changes do not transpose, rotate or reflect the field. The exported
image-reference origin and existing station indices/frame keys are retained.

Mean and instantaneous axial contour displays mask negative values and start
at zero. The saved data, statistics, profiles, distributions and quivers retain
all signed measurements. Transverse velocity, covariance and correlation
remain signed displays. This choice is a rendering policy, not sample rejection.
'''),nb.v4.new_code_cell('''from pathlib import Path
import sys,json
root=Path.cwd().resolve()
while not (root/'RAW_data').exists():
    if root==root.parent: raise RuntimeError('Workspace root not found')
    root=root.parent
sys.path.insert(0,str(root))
config_path=root/'PZT_Atomiser_Analysis/config/analysis.json'
config=json.loads(config_path.read_text())
assert config['velocity_scale']==0.290
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
