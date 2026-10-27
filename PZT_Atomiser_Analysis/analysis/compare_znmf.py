"""Compare raw export metadata without converting data or generating figures."""
import json
from pathlib import Path
import h5py
import numpy as np


def describe(path):
    with h5py.File(path, 'r') as f:
        attrs = {}
        datasets = 0
        def visit(name, obj):
            nonlocal datasets
            if obj.attrs:
                attrs[name] = {k:str(v) for k,v in obj.attrs.items()}
            if isinstance(obj,h5py.Dataset):
                datasets += 1
        visit('/',f)
        f.visititems(visit)
        keys = sorted(f['PIV FIELDS'],key=int)
        first = f['PIV FIELDS'][keys[0]]
        return {
            'path':str(path.resolve()),'bytes':path.stat().st_size,
            'root_groups':sorted(f.keys()),'dataset_count':datasets,
            'all_nonempty_attributes':attrs,
            'parameters':{k:np.asarray(v[()]).item() for k,v in f['PARAMETERS'].items() if not v.shape},
            'parameter_names':sorted(f['PARAMETERS'].keys()),
            'coordinates':{k:{'shape':list(f['PARAMETERS'][k].shape),
                              'dtype':str(f['PARAMETERS'][k].dtype),
                              'range':[float(f['PARAMETERS'][k][0]),float(f['PARAMETERS'][k][-1])],
                              'unique_steps':np.unique(np.diff(f['PARAMETERS'][k][:])).tolist()}
                           for k in ['x','y']},
            'fields':len(keys),'first_last_keys':[keys[0],keys[-1]],
            'first_field_schema':{k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in first.items()},
            'statistics_schema':{k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in f['STATISTICS'].items()},
            'NoSamples':int(f['STATISTICS/NoSamples'][()]),
        }


def main():
    root = Path(__file__).resolve().parents[2]
    znmf = root.parent/'ZNMF_Results'
    paths = {'PZT':root/'RAW_data/Batch4_Cam1.h5',
             'ZNMF_FOV1':znmf/'Raw data/First_location/Batch4_Cam1.h5',
             'ZNMF_FOV2':znmf/'Raw data/Second_location/Batch4_second_Cam1.h5',
             'ZNMF_FOV3':znmf/'Raw data/Third_location/Batch4_third_Cam1.h5'}
    result = {'files':{k:describe(p) for k,p in paths.items()}}
    pzt = result['files']['PZT']
    result['comparisons'] = {}
    for name,z in result['files'].items():
        if name == 'PZT':
            continue
        result['comparisons'][name] = {
            'same_root_groups':pzt['root_groups']==z['root_groups'],
            'same_parameter_names':pzt['parameter_names']==z['parameter_names'],
            'parameter_differences':{k:{'PZT':pzt['parameters'].get(k),'ZNMF':z['parameters'].get(k)}
                                     for k in set(pzt['parameters'])|set(z['parameters'])
                                     if pzt['parameters'].get(k)!=z['parameters'].get(k)},
            'same_field_dataset_names_and_dtypes':
                {k:v['dtype'] for k,v in pzt['first_field_schema'].items()} ==
                {k:v['dtype'] for k,v in z['first_field_schema'].items()},
            'same_statistics_names_and_dtypes':
                {k:v['dtype'] for k,v in pzt['statistics_schema'].items()} ==
                {k:v['dtype'] for k,v in z['statistics_schema'].items()},
            'same_y_coordinate_definition':pzt['coordinates']['y']==z['coordinates']['y'],
            'software_version_identity_verified':False,
            'run_export_settings_identity_verified':False,
            'flag_semantics_verified':False,
        }
    target = root/'PZT_Atomiser_Analysis/results/znmf_comparison.json'
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result['comparisons'],indent=2))
    print('Saved',target)


if __name__ == '__main__':
    main()
