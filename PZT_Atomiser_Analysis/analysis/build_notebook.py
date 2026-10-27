"""Compatibility entry point for the active calibrated notebook builder."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from PZT_Atomiser_Analysis.analysis.build_calibrated_notebook import main
if __name__ == '__main__':
    main()
