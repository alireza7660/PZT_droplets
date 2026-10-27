# Exact-vector comparison sheet

Configured conversion: V=(0.202) us + (0) vs and U=(0) us + (0.202) vs (m/s). Indices below are zero-based HDF5 [row,column]; frame is the exact group key, not a timestamp. Coordinates use the exported image-reference origin, not an outlet.

| Frame | [row,column] | Exported [x,y] | Axial x / transverse y (mm) | us / vs (assumed px/pair) | U / V (m/s) | Flag |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | [287,40] | [2591,512] | 2.585600 / 13.084550 | 0.49876404 / 3.8663788 | 0.78100851 / 0.10075034 | 1 |
| 1 | [159,40] | [2591,2560] | 12.928000 / 13.084550 | 1.2717133 / 3.8921127 | 0.78620677 / 0.25688608 | 1 |
| 1 | [32,40] | [2591,4592] | 23.189600 / 13.084550 | 0.5824216 / 1.4049326 | 0.28379639 / 0.11764916 | 1 |
| 879 | [287,40] | [2591,512] | 2.585600 / 13.084550 | 0.093841553 / 9.0510254 | 1.8283071 / 0.018955994 | 1 |
| 879 | [159,40] | [2591,2560] | 12.928000 / 13.084550 | 0.99848938 / 3.9863472 | 0.80524213 / 0.20169485 | 1 |
| 879 | [32,40] | [2591,4592] | 23.189600 / 13.084550 | -0.026292056 / 3.0809515 | 0.62235219 / -0.0053109953 | 1 |
| 1757 | [287,40] | [2591,512] | 2.585600 / 13.084550 | 0.66073608 / 3.6858826 | 0.74454828 / 0.13346869 | 1 |
| 1757 | [159,40] | [2591,2560] | 12.928000 / 13.084550 | -0.10408401 / 1.081234 | 0.21840926 / -0.021024971 | 1 |
| 1757 | [32,40] | [2591,4592] | 23.189600 / 13.084550 | -0.054236352 / 0.29632759 | 0.059858173 / -0.010955743 | 1 |

Full precision is retained in the adjacent vector_comparison.csv. Compare these exact nodes against a software numerical vector readout; if the software uses one-based array indices, add one to each row/column. Displayed coordinates may use a different reference, so check the exported coordinate fields as well.

A magnitude factor of 0.202 should apply only for pixel displacement storage. Check vertical sign using the same vector/component readout and the software axis legend, not the overall spray direction.