"""Compatibility entry point for the active calibrated notebook builder."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
if __name__ == '__main__':
    import json
    cfg=json.loads((Path(__file__).resolve().parents[1]/'config/analysis.json').read_text())
    if cfg.get('analysis_status')=='final':
        from PZT_Atomiser_Analysis.analysis.build_final_notebook import main
    else:
        from PZT_Atomiser_Analysis.analysis.build_calibrated_notebook import main
    main()
