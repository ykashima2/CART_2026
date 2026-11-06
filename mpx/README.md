# mpx: MPX batch 1

Molecular Pixelation (MPX) data of the first batch: eight samples from three donors, covering the CAR-T product (`t0`) and PBMCs after infusion (`t1`–`t3`). Sample information is in `data/sample_metadata.tsv`.

## Data

The `.pxl` files are not included in this repository. To replicate this analysis, place the PXL files for the samples listed in `data/sample_metadata.tsv` in the `data/` folder. The processing script matches each `sample` name in the metadata to a PXL file in `data/` (case-insensitive, searched recursively), so every file name must contain its sample name, e.g. `TKB_CAR-T_64` and `TKB_CAR_T_95` (note the different spelling in the sample names).

## Scripts

`scripts/01_processing_all_gating.qmd` processes the batch 1 samples and saves the processed object that `mpx_2/scripts/02_analysis_b1b2.qmd` uses for the manuscript figures. The other `.qmd` files in `scripts/` (t0/t1 and later-timepoint analyses, figure updates) are kept for the record and are not needed to reproduce the manuscript figures. The helper functions are in `scripts/R/`.

R packages are managed by the shared renv project in the root directory (`renv::restore()`, see the README there). `environment.yml` and `setup_conda_install_r_pkgs.sh` describe a separate conda setup that is not covered by `renv.lock`.
