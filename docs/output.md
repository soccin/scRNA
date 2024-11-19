# scRNA Pipeline Output

## Directory Layout

```
//seurat
├── stage1
│   ├── {PROJNO}_plt_01_QC.pdf
│   ├── {PROJNO}_plt_02_CellCycle.pdf
│   ├── {PROJNO}_plt_03_XenoStats.pdf
│   ├── {PROJNO}_plt_11_PostFilterQCTbls.pdf
│   ├── {PROJNO}_plt_11_PostFilterQCTbls.xlsx
│   ├── {PROJNO}_plt_12_PostIntegrateCC.pdf
│   └── {PROJNO}_plt_13_VariableFeatures.pdf
├── stage2
│   ├── {PROJNO}_plt_21_PCADimMetric.pdf
│   ├── {PROJNO}_plt_22_UMAP_20_.pdf
│   ├── {PROJNO}_plt_23_ClusterChart_20_.pdf
│   ├── {PROJNO}_plt_24_SampleBias_.pdf
│   └── {PROJNO}_plt_24_SampleBias_.xlsx
└── stage3
    ├── {PROJNO}_plt_31_ClusterMarkers_{RESOLUTION}_FDR_0.05_logFC_1.pdf
    |── {PROJNO}_plt_32_ClusterMarkersDot_{RESOLUTION}_FDR_0.05_logFC_1_.pdf
    ├── {PROJNO}_plt_33_ClusterHeatmap_{RESOLUTION}_FDR_0.05_logFC_1_.pdf
    └── {PROJNO}_plt_34_ClusterUMAP_{RESOLUTION}_FDR_0.05_logFC_1.pdf
```

