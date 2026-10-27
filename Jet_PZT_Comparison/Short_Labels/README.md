Short-axis-label version. Original files in the parent folder remain unchanged. Only the two axis labels differ in the copied plotting script; caption.txt defines s and the respective velocity components. Data exports and numerical verification are identical to the original.

Regenerate from this folder in PowerShell:

```powershell
python .\plot_comparison.py --znmf "..\..\..\ZNMF_Results\Final_Result\tables\FOV_1_2_3_centreline_and_width.csv" --pzt "..\..\PZT_Atomiser_Analysis\Final_Results\Batch4_Cam1\data\mean_v_vertical_profile.csv"
```

The 600 dpi PNG was visually inspected.
