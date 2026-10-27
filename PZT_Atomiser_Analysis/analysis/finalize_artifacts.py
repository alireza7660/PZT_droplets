"""Refresh notebook's embedded previews from current exported figure files."""
import base64
import json
from pathlib import Path
import shutil
import nbformat as nb


def main():
    root=Path(__file__).resolve().parents[2]
    cfg=json.loads((root/'PZT_Atomiser_Analysis/config/analysis.json').read_text())
    out=root/cfg['output']
    target=root/'PZT_Atomiser_Analysis/notebooks/Spatial_maps.ipynb'
    b=nb.read(target,as_version=4)
    b.cells[4].source=b.cells[4].source.replace("'frame_879_quiver'", "('frame_'+sorted(results[cid]['frames'],key=int)[len(results[cid]['frames'])//2]+'_quiver')")
    # Configuration is printed first and remains the central record if later edited.
    b.cells[0].source='''# PZT provisional calibrated analysis

Working assumptions, not final validated experimental results. User-confirmed
calibration is 5.05 micrometres/pixel and 25 microseconds. The source-component,
storage-unit, coordinate and sign assumptions are printed from the central
configuration below, and saved with every run. The initial provisional run
uses us/vs pixel displacements and retains their stored signs; those conventions
are not independently verified. All active results live in the configured
provisional_calibrated directory. The image-reference origin is not an outlet.
'''
    folder=out/'figures'/cfg['conditions'][0]['id']
    candidate=None;count=0
    for obj in b.cells[4].outputs:
        if obj.output_type=='stream':
            name=obj.text.strip()
            candidate=folder/(name+'.png')
        elif obj.output_type=='display_data' and 'image/png' in obj.data and candidate is not None and candidate.exists():
            obj.data['image/png']=base64.b64encode(candidate.read_bytes()).decode('ascii')
            candidate=None;count+=1
    for cell in b.cells:
        if cell.cell_type=='code':
            compile(cell.source,'notebook_cell','exec')
            assert cell.execution_count is not None
            assert not any(o.output_type=='error' for o in cell.outputs)
    nb.write(b,target)
    shutil.copy2(target,out/'notebooks'/target.name)
    print('Verified all code cells; synchronized embedded previews:',count)


if __name__=='__main__':main()
