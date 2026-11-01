"""Run the reusable comparison from this results folder or any working directory."""
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(root))
from PZT_Atomiser_Analysis.analysis.theory_comparison import run
print(run())
