"""Execute only measured plume alignment; preserve experimental cells and pause predictions."""
import asyncio
import hashlib
import json
import os
import shutil
import sys
import nbformat as nb
from nbclient import NotebookClient
from jupyter_client import AsyncKernelManager
from jupyter_client.kernelspec import KernelSpecManager
from .plume_alignment import ROOT,PROJECT,notebook_cells,json_out


def is_interaction(cell):
    return any(cell.metadata.get(k) for k in ('theory_comparison','combined_flow_prediction','prospective_studies','co_directed'))


def preserved(book):return [c for c in book.cells if not is_interaction(c) and not c.metadata.get('plume_alignment')]


def signature(cells):return hashlib.sha256(json.dumps(cells,sort_keys=True).encode()).hexdigest()


def main():
    cfg=json.loads((PROJECT/'config/plume_alignment.json').read_text());out=ROOT/cfg['output']
    target=PROJECT/'notebooks/PZT_final_analysis.ipynb';book=nb.read(target,as_version=4)
    original=signature(preserved(book));cells=notebook_cells()
    evidence=out/'validation';evidence.mkdir(parents=True,exist_ok=True)
    old_cells=[c for c in book.cells if is_interaction(c)]
    if old_cells:nb.write(nb.v4.new_notebook(cells=old_cells),evidence/'paused_interaction_notebook_section.ipynb')
    setup=nb.v4.new_code_cell('''from pathlib import Path
import sys,json
from IPython.display import display,Image
root=Path.cwd().resolve()
sys.path.insert(0,str(root))
output=root/'PZT_Atomiser_Analysis/Final_Results/Batch4_Cam1'
''')
    mini=nb.v4.new_notebook(cells=[setup,cells[1]])
    runtime=evidence/'jupyter_runtime';runtime.mkdir(exist_ok=True);os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3');spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'];manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(mini,timeout=600,resources={'metadata':{'path':str(ROOT)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager);client.execute()
    cells[1]=mini.cells[1];book.cells=preserved(book);book.cells[-1:-1]=cells
    assert signature(preserved(book))==original
    nb.write(book,target);shutil.copy2(target,PROJECT/'notebooks/Spatial_maps.ipynb')
    shutil.copy2(target,PROJECT/'Final_Results/Batch4_Cam1/methods/PZT_final_analysis.ipynb')
    json_out(evidence/'notebook_preservation.json',dict(original_noninteraction_cells_sha256=original,unchanged=True,
        executed_only='plume_alignment.run',interaction_predictions_run=False,previous_interaction_displays_preserved_separately=True))
    print('Executed measured plume alignment only; interaction predictions remain paused.')


if __name__=='__main__':main()
