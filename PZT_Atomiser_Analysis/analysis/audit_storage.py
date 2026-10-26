"""Read-only provenance audit; never infers units or flags from numeric values."""
import json
from pathlib import Path
import h5py


def audit(root):
    root = Path(root)
    source_files = [str(p.relative_to(root)) for p in root.rglob('*')
                    if p.is_file() and p.relative_to(root).parts[0]
                    not in {'.venv', 'PZT_Atomiser_Analysis', '.git'}]
    path = root/'RAW_data/Batch4_Cam1.h5'
    with h5py.File(path, 'r') as f:
        objects, attributes, types = [], {}, set()
        def visit(name, obj):
            objects.append(name)
            if len(obj.attrs):
                attributes[name] = {k:str(v) for k,v in obj.attrs.items()}
            if isinstance(obj,h5py.Dataset):
                types.add(str(obj.dtype))
        f.visititems(visit)
        result = {
            'original_workspace_files_excluding_generated_analysis_and_environment':sorted(source_files),
            'hdf5_path':str(path.relative_to(root)),
            'root_groups':list(f), 'root_attributes':dict(f.attrs),
            'userblock_size':f.userblock_size,'objects_visited':len(objects),
            'nonempty_object_attributes':attributes, 'dataset_dtypes':sorted(types),
            'parameters':{k:f['PARAMETERS'][k][()].item() for k in f['PARAMETERS']
                          if f['PARAMETERS'][k].shape == ()},
            'coordinates':{k:{'first':float(f['PARAMETERS/'+k][0]),
                               'last':float(f['PARAMETERS/'+k][-1]),
                               'first_step':float(f['PARAMETERS/'+k][1]-f['PARAMETERS/'+k][0])}
                           for k in ['x','y']},
            'exporter_source_found':False,
            'storage_units_verified':False,'component_basis_verified':False,
            'prior_transformations_verified':False,'flag_definitions_verified':False,
        }
    dest = root/'PZT_Atomiser_Analysis/results/storage_evidence.json'
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    audit(Path(__file__).resolve().parents[2])
