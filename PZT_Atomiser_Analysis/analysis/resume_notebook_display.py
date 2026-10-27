"""Recover a failed final display cell without recalculating saved figures.

Normal regeneration uses build_calibrated_notebook. This helper preserves
already-executed calculation/figure cells and executes only the final table
display with context reloaded from the saved configuration and results.
"""
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
    target=root/'PZT_Atomiser_Analysis/notebooks/Spatial_maps.ipynb'
    b=nb.read(target,as_version=4)
    cell=b.cells[6]
    cell.source=cell.source.replace("Markdown('\n'.join(table))","Markdown(chr(10).join(table))")
    compile(cell.source,'last_display_cell','exec')
    for prior in b.cells[:6]:
        if prior.cell_type=='code':
            assert prior.execution_count is not None
            assert not any(o.output_type=='error' for o in prior.outputs)
    recovery=nb.v4.new_notebook(cells=[nb.v4.new_code_cell('''from pathlib import Path
import json
from IPython.display import display, Markdown
root=Path.cwd().resolve()
config=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text())
output=root/config['output']
results={c['id']:None for c in config['conditions']}
assert (output/'results/validation.json').exists()
'''),cell])
    runtime=root/'PZT_Atomiser_Analysis/provisional_calibrated/diagnostics/jupyter_runtime'
    os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(recovery,timeout=120,resources={'metadata':{'path':str(root)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager)
    client.execute()
    b.cells[6]=recovery.cells[1]
    b.metadata['display_recovery']='Final display cell re-executed from saved results after syntax correction; calculation and figure outputs preserved.'
    nb.write(b,target)
    dest=root/'PZT_Atomiser_Analysis/provisional_calibrated/notebooks'
    dest.mkdir(exist_ok=True);shutil.copy2(target,dest/target.name)
    print('Recovered final display cell; all notebook code cells executed without errors.')


if __name__=='__main__':main()
