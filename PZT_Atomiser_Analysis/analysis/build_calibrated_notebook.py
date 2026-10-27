"""Build and execute the notebook that regenerates the provisional result set."""
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
    out=root/'PZT_Atomiser_Analysis'
    runtime=out/'provisional_calibrated/diagnostics/jupyter_runtime'
    runtime.mkdir(parents=True,exist_ok=True)
    os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    b=nb.v4.new_notebook()
    b.metadata.kernelspec={'display_name':'Python (PZT workspace)','language':'python','name':'python3'}
    b.cells=[nb.v4.new_markdown_cell('''# PZT provisional calibrated analysis

**Working assumptions, not final validated experimental results.** Calibration
is user-confirmed: 5.05 micrometres/pixel and 25 microseconds. Source us/vs are
provisionally horizontal/vertical pixel displacements. Stored vertical sign
is retained, so U=+0.202 vs and V=+0.202 us. Exporter units, signs, upstream
replacement/smoothing and flag meanings remain unverified. No new filtering,
gap filling or extrapolation is applied. The image-reference origin is not
an established outlet. All result files live in provisional_calibrated.
'''),nb.v4.new_code_cell('''from pathlib import Path
import sys, json
root=Path.cwd().resolve()
while not (root/'RAW_data').exists():
    if root==root.parent: raise RuntimeError('Workspace root not found')
    root=root.parent
sys.path.insert(0,str(root))
config_path=root/'PZT_Atomiser_Analysis/config/analysis.json'
config=json.loads(config_path.read_text())
print(json.dumps(config,indent=2))
'''),nb.v4.new_code_cell('''from PZT_Atomiser_Analysis.analysis.pipeline import run
results=run(config_path)
output=root/config['output']
'''),nb.v4.new_markdown_cell('''## Methods and population

The reusable pipeline scales instantaneous components once, then computes
means, RMS, speed and signed covariance from the same finite-pair population.
Denominator N is used. No flag is rejected. Flag-specific means/RMS are
selection-sensitivity diagnostics, not physically validated populations.
Full-range colour limits are shared across matching quantities and flag groups.
Profiles use native rows. Each station PDF is a distinct native-node empirical
distribution; histogram PDFs and exact empirical CDFs preserve sample support.
No independent-frame uncertainty or turbulence interpretation is claimed.
'''),nb.v4.new_code_cell('''from IPython.display import Image, Markdown, display
display(Markdown((output/'METHODS.md').read_text(encoding='utf-8')))
for cid in results:
    folder=output/'figures'/cid
    for name in ['mean_U','mean_V','magnitude_mean_vector','mean_instantaneous_speed',
                 'rms_U','rms_V','covariance','correlation','sampling_stations',
                 'profile_mean_U','profile_rms_U','U_pdf','V_pdf','speed_pdf',
                 'U_ecdf',('frame_'+sorted(results[cid]['frames'],key=int)[len(results[cid]['frames'])//2]+'_quiver')]:
        print(name)
        display(Image(filename=str(folder/(name+'.png')),width=600))
'''),nb.v4.new_markdown_cell('''## Exact vectors for software comparison

All row/column indices are zero-based HDF5 array indices. Compare the same
numeric frame key and node in the processing software, including the displayed
coordinate reference and component units. These samples can test the 0.202
factor and retained vertical sign without inferring either from mean direction.
'''),nb.v4.new_code_cell('''import csv
for cid in results:
    with (output/'results'/f'{cid}_vector_comparison.csv').open() as f:
        rows=list(csv.DictReader(f))
    columns=list(rows[0])
    table=[' | '.join(columns), ' | '.join(['---']*len(columns))]
    table += [' | '.join(str(row[c]) for c in columns) for row in rows]
    display(Markdown(chr(10).join(table)))
    with (output/'results'/f'{cid}_flag_sensitivity_summary.csv').open() as f:
        print(f.read())
    display(Markdown((output/'FINDINGS.md').read_text(encoding='utf-8')))
    print((output/'results/validation.json').read_text())
'''),nb.v4.new_markdown_cell('''## Regenerating after software verification

Edit components, component_mapping, coordinate_mapping, transform, scales,
origin and working_assumptions in config/analysis.json. The plot basis is
[transverse V, axial U]. A common transform is applied to both positions and
velocities. An exporter-only component correction is distinct from a physical
basis change; record the evidence rather than reflecting coordinates twice.
Rerun this notebook to regenerate every result, including covariance, arrows,
profiles, station samples and the vector comparison sheet. No derivatives are
computed in this run. Configuration and methods are saved beside the results.
''')]
    target=out/'notebooks/Spatial_maps.ipynb';target.parent.mkdir(exist_ok=True)
    nb.write(b,target)
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(b,timeout=1200,resources={'metadata':{'path':str(root)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager)
    try:client.execute()
    finally:nb.write(b,target)
    dest=out/'provisional_calibrated/notebooks';dest.mkdir(exist_ok=True)
    shutil.copy2(target,dest/target.name)
    print('Executed notebook:',target)


if __name__=='__main__':main()
