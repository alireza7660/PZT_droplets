"""Regenerate the input-gated combined workflow and only its notebook section."""
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
from .combined_flow_prediction import ROOT,PROJECT,notebook_cells,write_json


def unrelated_cells(book):
    return [c for c in book.cells if not c.metadata.get('theory_comparison') and not c.metadata.get('combined_flow_prediction') and not c.metadata.get('prospective_studies')]


def signature(cells):return hashlib.sha256(json.dumps(cells,sort_keys=True).encode()).hexdigest()


def main():
    cfg=json.loads((PROJECT/'config/combined_flow_prediction.json').read_text());out=ROOT/cfg['output']
    target=PROJECT/'notebooks/PZT_final_analysis.ipynb';book=nb.read(target,as_version=4)
    original=signature(unrelated_cells(book));cells=notebook_cells()
    setup=nb.v4.new_code_cell('''from pathlib import Path
import sys,json
from IPython.display import display,Image
root=Path.cwd().resolve()
sys.path.insert(0,str(root))
output=root/'PZT_Atomiser_Analysis/Final_Results/Batch4_Cam1'
''')
    mini=nb.v4.new_notebook(cells=[setup,cells[1],cells[3]])
    runtime=out/'validation/jupyter_runtime';runtime.mkdir(parents=True,exist_ok=True)
    os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'];manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(mini,timeout=600,resources={'metadata':{'path':str(ROOT)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager);client.execute()
    cells[1]=mini.cells[1];cells[3]=mini.cells[2];book.cells=unrelated_cells(book);book.cells[-1:-1]=cells
    assert signature(unrelated_cells(book))==original
    nb.write(book,target);shutil.copy2(target,PROJECT/'notebooks/Spatial_maps.ipynb')
    shutil.copy2(target,ROOT/'PZT_Atomiser_Analysis/Final_Results/Batch4_Cam1/methods/PZT_final_analysis.ipynb')
    write_json(out/'validation/notebook_preservation.json',dict(nonprediction_cells_sha256=original,experimental_cells_preserved=True))
    from .validate_combined_flow import main as validate
    validate()
    if cfg.get('design_study',{}).get('enabled'):
        from .validate_prospective_studies import main as validate_studies
        validate_studies()
    print('Updated and executed only the combined-flow workflow notebook section.')


if __name__=='__main__':main()
