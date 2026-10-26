import json
from pathlib import Path
import h5py
import numpy as np


def inspect(path):
    inventory = {}
    def clean(value):
        if isinstance(value, bytes):
            return value.decode('utf-8', errors='replace')
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        return value
    with h5py.File(path, 'r') as f:
        inventory['/'] = {'attributes': {k: clean(v) for k,v in f.attrs.items()}}
        def visit(name, obj):
            item = {'attributes': {k: clean(v) for k,v in obj.attrs.items()}}
            if isinstance(obj, h5py.Dataset):
                item.update(shape=obj.shape, dtype=str(obj.dtype))
                if obj.size < 30:
                    item['values'] = clean(obj[()])
            inventory[name] = item
        f.visititems(visit)
    return inventory


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    result = inspect(root/'RAW_data/Batch4_Cam1.h5')
    target = root/'PZT_Atomiser_Analysis/results/input_inventory.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, default=str), encoding='utf-8')
    print(f'Saved {len(result)} inventory entries to {target}')
