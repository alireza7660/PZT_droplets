"""Execute only the added profile cell; preserve existing notebook calculations."""
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
from .vertical_profile import notebook_cells


def main():
    root=Path(__file__).resolve().parents[2];project=root/'PZT_Atomiser_Analysis'
    cfg=json.loads((project/'config/analysis.json').read_text());out=root/cfg['output']
    target=project/'notebooks/PZT_final_analysis.ipynb';book=nb.read(target,as_version=4)
    cells=notebook_cells()
    setup=nb.v4.new_code_cell('''from pathlib import Path
import sys,json
from IPython.display import display,Image
root=Path.cwd().resolve()
sys.path.insert(0,str(root))
config=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text())
output=root/config['output']
''')
    mini=nb.v4.new_notebook(cells=[setup,cells[1]])
    runtime=out/'diagnostics/jupyter_runtime';runtime.mkdir(exist_ok=True)
    os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'];manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(mini,timeout=120,resources={'metadata':{'path':str(root)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager);client.execute()
    # Idempotently replace only previously added profile cells, if present.
    book.cells=[c for c in book.cells if not c.get('metadata',{}).get('fixed_x_vertical_profile')]
    cells[1]=mini.cells[1]
    for c in cells:c.metadata['fixed_x_vertical_profile']=True
    book.cells[-1:-1]=cells
    # Synchronize displayed summary text, whose only changes are the new figure
    # count and additional selection metadata; existing figure payloads are retained.
    for cell in book.cells:
        if cell.cell_type!='code':continue
        for obj in cell.outputs:
            if obj.output_type=='stream' and 'Measured extrema and display limits are reported separately:' in obj.text:
                obj.text='Measured extrema and display limits are reported separately:\n'+(out/'methods/summary.json').read_text()+'\n'+(out/'methods/validation.json').read_text()+'\n'
            elif obj.output_type=='display_data' and obj.get('data',{}).get('text/markdown','').startswith('# PZT atomiser methods'):
                obj.data['text/markdown']=(out/'methods/METHODS.md').read_text(encoding='utf-8')
    nb.write(book,target);shutil.copy2(target,project/'notebooks/Spatial_maps.ipynb');shutil.copy2(target,out/'methods'/target.name)
    print('Added and executed the fixed-x profile cell; existing calculation cells and figures retained.')


if __name__=='__main__':main()
