From: (https://btep.ccr.cancer.gov/question/single_cell_rna_seq/can-you-give-some-suggestions-cell-identity-annoation-and-recommend-packages-good-at-cell-identity-annotation-what-is-the-most-used-practice-manually-or-auto/) 

Seurat also provides an additional option for cell type identification with its AddModuleScore function. This approach is implemented by providing gene sets characteristic of different cell types and letting Seurat compute a score for each cell type for all cells in the data. Using this approach, the highest scoring cell type (per cell) is assigned as the cell type. Seurat actually uses the very same AddModuleScore function for mapping cells to different cell cycle phases by utilizing the canonical markers of G1, S, G2/M phases.

