"""Refresh signed contour colour-bar ticks and matching embedded previews."""
import base64
import json
from pathlib import Path
import shutil
import numpy as np
import nbformat as nb
from .plotting import contour_map
from .final_analysis import AXES,LABELS,UNIT


def main():
    root=Path(__file__).resolve().parents[1]
    cfg=json.loads((root/'config/analysis.json').read_text());out=root.parent/cfg['output']
    manifest=json.loads((out/'methods/figure_manifest.json').read_text());replacements={}
    for entry in manifest:
        key=entry['id']
        signed=key in ['mean_u','covariance','correlation'] or (key.startswith('frame_') and key.endswith('_u')) or (key.startswith('flag_') and key.endswith('_mean_u'))
        if not signed:continue
        png=out/entry['files']['png'];old=base64.b64encode(png.read_bytes()).decode('ascii')
        data=out/entry['numerical_source']
        with np.load(data) as z:
            x,y,values,lim=z['x'].copy(),z['y'].copy(),z['values'].copy(),z['color_limits'].copy()
            levels=len(z['levels'])
        label=rf'$u$ ({UNIT})' if key.startswith('frame_') else LABELS['mean_u' if key.startswith('flag_') else key]
        contour_map(x,y,values,label=label,limits=lim,output=png.with_suffix(''),**AXES,signed=True,levels=levels,data_output=data)
        replacements[old]=base64.b64encode(png.read_bytes()).decode('ascii')
    target=root/'notebooks/PZT_final_analysis.ipynb';b=nb.read(target,as_version=4);count=0
    for cell in b.cells:
        if cell.cell_type!='code':continue
        assert cell.execution_count is not None
        assert not any(o.output_type=='error' for o in cell.outputs)
        for obj in cell.outputs:
            if obj.output_type in ['display_data','execute_result'] and obj.get('data',{}).get('image/png') in replacements:
                obj.data['image/png']=replacements[obj.data['image/png']];count+=1
    nb.write(b,target);shutil.copy2(target,root/'notebooks/Spatial_maps.ipynb');shutil.copy2(target,out/'methods'/target.name)
    print('Refreshed',len(replacements),'signed contour sets and',count,'embedded previews.')


if __name__=='__main__':main()
