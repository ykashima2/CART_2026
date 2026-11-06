# mpx_2: MPX batches 1 + 2

Molecular Pixelation (MPX) data of the second batch (donors D–H) and the combined analysis of both batches (t0 and t1, 14 samples, 8 donors). This project generates the MPX manuscript figures. Sample information for batch 2 is in `data/sample_metadata.tsv`.

## Data

The `.pxl` files are not included in this repository. To replicate this analysis, place the PXL files listed in the `file_path` column of `data/sample_metadata.tsv` in this project's `data/` folder (`TKB_CART_<nr>.layout.dataset.pxl`), and the batch 1 PXL files listed in `mpx/data/sample_metadata.tsv` in `mpx/data/`. The batch 1 files are also read by `scripts/01_processing_t0.qmd` and by the cell graph panels of `scripts/02_analysis_b1b2.qmd`.

## Scripts and run order

1. `mpx/scripts/01_processing_all_gating.qmd` (batch 1)
2. `scripts/01_processing_all_gating.qmd` (batch 2)
3. `scripts/01_processing_t0.qmd` (CAR-T product, both batches)
4. `scripts/02_analysis_b1b2.qmd` (combined analysis; panels flagged `Manuscript:` are the manuscript figures, written to `results/b1b2_t0t1/revised_figs/`)

The helper functions are in `scripts/R/`. The other `.qmd` files in `scripts/` are kept for the record and are not needed to reproduce the manuscript figures. The PNA manuscript figures (`pna_b1b2/modules/08_manuscript_figures.qmd`) use results from step 4, so run it first.

R packages are managed by the shared renv project in the root directory (`renv::restore()`, see the README there).
