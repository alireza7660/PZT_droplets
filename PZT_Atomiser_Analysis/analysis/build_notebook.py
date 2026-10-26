"""Create and execute the notebook using this environment's Python kernel."""
from pathlib import Path
import os
import sys
import nbformat as nb
from nbclient import NotebookClient
from jupyter_client.kernelspec import KernelSpecManager

root = Path(__file__).resolve().parents[2]
out = root/'PZT_Atomiser_Analysis'
runtime = out/'diagnostics/jupyter_runtime'
runtime.mkdir(parents=True,exist_ok=True)
os.environ['JUPYTER_RUNTIME_DIR'] = str(runtime)
b = nb.v4.new_notebook()
b.metadata.kernelspec = {'display_name':'Python (PZT workspace)','language':'python','name':'python3'}
b.cells = [nb.v4.new_markdown_cell('''# PZT spatial fields

This notebook calculates statistics and generates actual PDF, SVG and 600-dpi
PNG figures through reusable functions. Spatial calibration (5.05 um/pixel)
and pulse separation (25 us) are user-confirmed. Stored-array units, axial
orientation, flags and the processing/export convention remain unverified.
Plots are therefore explicitly uncalibrated native-component diagnostics.
No filtering, gap filling or extrapolation is applied.
'''),nb.v4.new_code_cell('''from pathlib import Path
import sys
root = Path.cwd().resolve()
while not (root / 'RAW_data').exists():
    if root == root.parent:
        raise RuntimeError('Workspace root not found')
    root = root.parent
sys.path.insert(0, str(root))
from PZT_Atomiser_Analysis.analysis.pipeline import run
results = run(root / 'PZT_Atomiser_Analysis/config/analysis.json')
'''),nb.v4.new_markdown_cell('''## Definitions and validation

Both components use paired finite samples; denominator N is used for RMS and
covariance. No flag meanings are inferred. All current u/v pairs are finite,
so finite support does not establish quality acceptance. Saved U/V match means
of us/vs, not u/v. The current analysis recomputes all moments from u/v.
Independent point checks, variance and correlation bounds, and the mean-speed
inequality are checked before plotting. Contours mask any grid cell touching
a missing node. Full-range limits are shared across configured conditions;
instantaneous limits are shared across selected frames. The first, middle and
last numeric field keys are selected without implying phase or typicality.
'''),nb.v4.new_code_cell('''from IPython.display import display, Image
import json
out = root / 'PZT_Atomiser_Analysis'
for cid, r in results.items():
    print(cid, json.dumps(r['metadata'], indent=2))
    folder = out / 'diagnostics/uncalibrated' / cid
    for quantity in ['mean_u', 'mean_v', 'magnitude_mean_vector',
                     'mean_instantaneous_speed', 'rms_u', 'rms_v',
                     'covariance', 'correlation', 'finite_pair_fraction',
                     'flag_0_fraction']:
        print(quantity)
        display(Image(filename=str(folder / (quantity + '.png')), width=650))
'''),nb.v4.new_markdown_cell('''## Reuse and remaining definitions

Edit config/analysis.json for another condition and rerun the notebook.
Coordinate and velocity scale factors must follow the verified exporter
convention. The factor 0.202 applies only if u/v are pixel displacements per
pair. Physical axial/transverse mapping requires original image orientation.
Flag selection requires a documented definition source. See ../README.md.
''')]
target = out/'notebooks/Spatial_maps.ipynb'
target.parent.mkdir(exist_ok=True)
nb.write(b,target)
client = NotebookClient(b, timeout=900, resources={'metadata':{'path':str(root)}})
# Specify the workspace interpreter without installing a global kernelspec.
ksm = KernelSpecManager()
spec = ksm.get_kernel_spec('python3')
spec.argv = [sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
ksm.get_kernel_spec = lambda name: spec
from jupyter_client import AsyncKernelManager
client.km = AsyncKernelManager(kernel_name='python3',kernel_spec_manager=ksm)
try:
    client.execute()
finally:
    nb.write(b,target)
print('Executed notebook:',target)
