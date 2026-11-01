"""Execute only the new theory section; retain prior experimental notebook outputs."""
import asyncio
import os
from pathlib import Path
import shutil
import sys
import nbformat as nb
from nbclient import NotebookClient
from jupyter_client import AsyncKernelManager
from jupyter_client.kernelspec import KernelSpecManager
from .theory_comparison import notebook_cells,ROOT,PROJECT,FINAL


def main():
    target=PROJECT/'notebooks/PZT_final_analysis.ipynb'
    book=nb.read(target,as_version=4);cells=notebook_cells()
    setup=nb.v4.new_code_cell('''from pathlib import Path
import sys,json
from IPython.display import display,Image
root=Path.cwd().resolve()
sys.path.insert(0,str(root))
output=root/'PZT_Atomiser_Analysis/Final_Results/Batch4_Cam1'
''')
    mini=nb.v4.new_notebook(cells=[setup,cells[1]])
    runtime=FINAL/'Theory_Comparison/jupyter_runtime';runtime.mkdir(exist_ok=True)
    os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
    if sys.platform=='win32':asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    manager=KernelSpecManager();spec=manager.get_kernel_spec('python3')
    spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    manager.get_kernel_spec=lambda name:spec
    client=NotebookClient(mini,timeout=300,resources={'metadata':{'path':str(ROOT)}})
    client.km=AsyncKernelManager(kernel_name='python3',kernel_spec_manager=manager)
    client.execute()
    cells[1]=mini.cells[1]
    book.cells=[c for c in book.cells if not c.metadata.get('theory_comparison')]
    book.cells[-1:-1]=cells
    nb.write(book,target)
    shutil.copy2(target,PROJECT/'notebooks/Spatial_maps.ipynb')
    shutil.copy2(target,FINAL/'methods/PZT_final_analysis.ipynb')
    print('Executed added theory section; existing experimental cells retained.')


if __name__=='__main__':main()
