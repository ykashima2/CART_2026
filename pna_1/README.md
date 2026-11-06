<div align="center">

# Proxiome Analysis Template

**A complete, reproducible workflow for single-cell Proximity Network Analysis**

[![R](https://img.shields.io/badge/R-4.4.3-blue?logo=r)](https://www.r-project.org/)
[![Quarto](https://img.shields.io/badge/Quarto-Document-purple)](https://quarto.org/)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-green.svg)](LICENSE.md)

*Analyze PNA data from the [nf-core/pixelator](https://nf-co.re/pixelator) pipeline with quality control, clustering, annotation, and publication-ready visualizations.*

[**Get Started →**](#quick-start) · [Tutorials](https://software.pixelgen.com/pna-analysis/introduction/) · [pixelatorR](https://github.com/PixelgenTechnologies/pixelatorR)

</div>

---

<br>

## Contents

- [Sample Information](#sample-information)
- [Quick Start](#quick-start)
- [First-Time R Setup](#first-time-r-setup)
- [Prepare Your Data](#prepare-your-data)
- [Run the Analysis](#run-the-analysis)
- [Performed analysis description](#performed-analysis-description)
- [Troubleshooting](#troubleshooting)
- [Maintainers](#maintainers)

---

<br>

## Sample Information

This analysis covers **8 PNA samples** from **5 donors**, all carrying the **tisa-cel** CAR construct and all from responding donors. Samples were collected either from the CAR-T product itself (`t0`, day 0) or from PBMCs drawn after infusion (`t1`, 7–11 days post-infusion).

| sample_id | sample_alt_name | sample_alias | dpi | donor | cell_viability | condition | file_path | timepoint | construct | status |
|-----------|-----------------|--------------|:---:|-------|:--------------:|-----------|-----------|:---------:|-----------|--------|
| TKB_CAR-T_10 | TKB_CAR-T_10 | donorD_t1 | 7 | donorD | 76 | PBMC | TKB_CAR-T_10.layout.pxl | t1 | tisa-cel | responding |
| TKB_CAR-T_2 | TKB_CAR-T_2 | donorE_t1 | 7 | donorE | 59 | PBMC | TKB_CAR-T_2.layout.pxl | t1 | tisa-cel | responding |
| TKB_CAR-T_50 | TKB_CAR-T_50 | donorF_t0 | 0 | donorF | 9 | CAR-T_PRODUCT | TKB_CAR-T_50.layout.pxl | t0 | tisa-cel | responding |
| TKB_CAR-T_55 | TKB_CAR-T_55 | donorF_t1 | 11 | donorF | 95 | PBMC | TKB_CAR-T_55.layout.pxl | t1 | tisa-cel | responding |
| TKB_CAR-T_194 | TKB_CAR-T_194 | donorG_t0 | 0 | donorG | 29 | CAR-T_PRODUCT | TKB_CAR-T_194.layout.pxl | t0 | tisa-cel | responding |
| TKB_CAR-T_198 | TKB_CAR-T_198 | donorG_t1 | 7 | donorG | 87 | PBMC | TKB_CAR-T_198.layout.pxl | t1 | tisa-cel | responding |
| TKB_CAR-T_175 | TKB_CAR-T_175 | donorH_t0 | 0 | donorH | 37 | CAR-T_PRODUCT | TKB_CAR-T_175.layout.pxl | t0 | tisa-cel | responding |
| TKB_CAR-T_181 | TKB_CAR-T_181 | donorH_t1 | 10 | donorH | 96 | PBMC | TKB_CAR-T_181.layout.pxl | t1 | tisa-cel | responding |

**At a glance**

| | |
|---|---|
| Samples | 8 |
| Donors | 5 (donorD–donorH) |
| Conditions | `CAR-T_PRODUCT` (3) · `PBMC` (5) |
| Timepoints | `t0` (3) · `t1` (5) |
| Days post-infusion (`dpi`) | 0 (3) · 7 (3) · 10 (1) · 11 (1) |
| Paired donors (t0 + t1) | 3 of 5 (donorF, donorG, donorH) |
| Construct | tisa-cel (8/8) |
| Status | responding (8/8) |
| Cell viability | 9–96 % (median 67.5 %) |

> 📝 Donors D and E contribute a `t1` sample only, so 3 of 5 donors are paired across timepoints. The CAR-T product samples have markedly lower viability (9–37 %) than the post-infusion PBMCs (59–96 %). `sample_alt_name` duplicates `sample_id` in this dataset.

> 📁 **PXL files.** The `.pxl` files are not included in this repository. To replicate this analysis, place the PXL files listed in `data/metadata.csv` (`<sample_id>.layout.pxl`) in the `data/` folder of this project.

> ⚠️ This table mirrors `data/metadata.csv`. If you edit the CSV (samples, order, values), update this section to match. Each `sample_id` corresponds to `<sample_id>.layout.pxl` in `data/`.

---

<br>

## Quick Start

This folder is one of several analysis projects in the same repository, so there is nothing to download or clone separately. Continue to [First-Time R Setup](#first-time-r-setup).

> 💡 **Tip:** Always open your project by double-clicking `proxiome_analysis_template.Rproj` (found directly in this folder) - this sets up file paths automatically.

---

<br>

## First-Time R Setup

### 1. Install R and RStudio

You need to download and install **both** programs:

| Software | Description | Download |
|----------|-------------|----------|
| **R** (4.4.3) | The programming language | [Download from CRAN](https://cran.r-project.org/) |
| **RStudio** | The editor/interface for R | [Download from Posit](https://posit.co/download/rstudio-desktop/) |

> 🆕 **New to R?** Read the [Getting Started with R](https://rstudio-education.github.io/hopr/starting.html) guide first.

<br>

### 2. Install Compiler Tools

Some R packages require compilation. Install the appropriate tools for your operating system:

<br>

<details>
<summary><strong>Windows</strong></summary>

Also install [Rtools](https://cran.r-project.org/bin/windows/Rtools/). Choose the version matching your R installation (required for some packages).

</details>

<details>
<summary><strong>macOS</strong></summary>

Press Command + Spacebar to open Spotlight search, type "Terminal," and hit Enter to open the terminal.

Install Xcode command line tools by running this in the terminal: 

```bash
xcode-select --install
```

> ⚠️ **Note:** This command will open a separate installation window. Follow the prompts in that window to complete the installation before continuing.

</details>

<details>
<summary><strong>Linux</strong></summary>

After installing R, install system libraries: 

```bash
sudo apt install cmake libglpk-dev libhdf5-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev
```

</details>

<br>

### 3. Restore the R packages

Package versions are pinned in the shared renv project at the root directory (`renv.lock`: R 4.4.3, with `pixelatorR` installed from GitHub). Open the R project file in RStudio by double-clicking `proxiome_analysis_template.Rproj` (located directly in this folder). If renv asks to restore the library, or `renv::status()` reports differences, run the following in the R console:

```r
renv::restore()
```

For details, see the README in the root directory.

---

<br>

## Prepare Your Data

When you have everything installed, add your `.pxl` data files and a `metadata.csv` to the `data/` folder as shown below.

### Project structure

```
pna_1/
├── data/
│   ├── metadata.csv                  ← Sample info
│   └── *.layout.pxl                  ← PXL files
├── modules/                          ← Analysis notebooks, all samples (run in order)
│   ├── common_setup.R                ← Shared setup (sourced automatically)
│   ├── 01_quality_control.qmd        ← Start here (load data + QC)
│   ├── 02_clustering_annotation.qmd
│   ├── 03_abundance.qmd
│   ├── 04_proximity.qmd
│   └── 05_statistical_testing.qmd    ← Optional
├── modules_t0/                       ← Same notebooks for the t0 product-only sub-analysis
├── results/                          ← Created automatically when running the modules
├── results_t0/                       ← Output of modules_t0/
├── LICENSE.md
├── proxiome_analysis_template.Rproj  ← Double-click to open the project in RStudio
├── proxiome_analysis_template.qmd    ← Full re-render of the modules (optional)
└── proxiome_analysis_template_t0.qmd ← Full re-render of the t0 sub-analysis (optional)
```

### Create metadata.csv

The `metadata.csv` file describes **your** samples and links them to the corresponding PXL files. 

| Column | Required | Description | Example |
|--------|:--------:|-------------|---------|
| `sample_id` | ✓ | Unique identifier | `S1` |
| `sample_alias` | ✓ | Descriptive name | `S1_donor1` |
| `condition` | ✓ | Experimental group | `donor1` |
| `file_path` | ✓ | PXL filename (not full path) | `Sample1.layout.pxl` |
| *other* | | Any extra columns (`donor`, `batch`, etc.) | - |

> 💡 Extra columns are available as variables in the template for grouping or covariates.

**Example file:**

The release zip ships `data/metadata.csv` with only the header row. Add one row per sample. In this repository, `data/metadata.csv` is filled for donors 1 and 3 from the public [Resting PBMCs from 16 healthy donors](https://software.pixelgen.com/datasets/PBMC-16-healthy-donors-v2.0-proxiome-immuno-155/) dataset (both samples are resting; `condition` is set to donor identity so module `05` can contrast them):

```csv
sample_id,sample_alias,condition,file_path,donor,age,sex
S1,S1_donor1,donor1,PBMCs_resting_S01_HG2.layout.pxl,donor1,58,Female
S2,S2_donor3,donor3,PBMCs_resting_S01_HG6.layout.pxl,donor3,21,Female
```

> 📝 Row order controls plot order. Both `,` and `;` separators work. Every `file_path` must match a file in `data/`.

---

<br>

## Run the Analysis

The analysis lives in numbered notebooks under `modules/` — **Quarto documents** (QMD files) that combine text explanations with runnable R code. [PNA tutorials](https://software.pixelgen.com/pna-analysis/introduction/) provide additional background.

- **QMD file**: A document that mixes formatted text with executable R code. You run the code section by section.
- **Code chunk**: A gray block containing R code. Each chunk has a green **▶ (play) button** in the top-right corner.

### How to run the analysis (first time)

1. **Open** `proxiome_analysis_template.Rproj`
2. **Open** `modules/01_quality_control.qmd` in RStudio
3. **Run chunks** top to bottom through module `01`, then continue with `02`, `03`, and `04`
4. **Optional:** Open `modules/05_statistical_testing.qmd` for differential testing

In this project modules `01`–`05` were run on all samples. The product-only sub-analysis of the `t0` samples uses the notebooks in `modules_t0/` (modules `01`–`04`), and `proxiome_analysis_template_t0.qmd` re-renders it.

> 💡 **Tip:** Each module starts with a short intro. Search for **Decision checkpoint** to find values you need to edit.

### Iterate on downstream modules

After module `02` finishes, you can open `03`, `04`, or `05` in a **new R session** and run the first chunk — packages, palettes, and saved data load automatically.

| Checkpoint file | Created by | Used by |
|-----------------|------------|---------|
| `pg_data_merged.rds` | module `01` (filtered cells) | module `02` |
| `selected_markers.rds` | module `01` | module `02` |
| `annotated_seurat_object.rds` | module `02` | modules `03`–`05` |
| `markers_to_test.rds` | module `02` | modules `04`–`05` |

All checkpoints are saved under `results/checkpoint_data/`.

### Full re-render (optional)

To regenerate one HTML report of the entire workflow, install [Quarto](https://quarto.org/) and render `proxiome_analysis_template.qmd`. For day-to-day analysis, use the individual modules instead.

### ⚠️ Where you need to make choices

| Module | What to do |
|--------|------------|
| **01 — sample table** | Check `metadata.csv` paths and sample IDs match your data |
| **01 — QC thresholds** | Set `n_umi_min_threshold`, `n_umi_max_threshold`, `isotype_percent_threshold` after reviewing QC plots; re-run section 1.3.5 onward |
| **01 — marker selection** | Review `selected_markers` after running the selection code |
| **02 — manual annotation** | **Critical step!** Edit `cluster_cell_annotation` to assign cell type names to clusters |
| **05 — reference condition** (optional) | Set `reference_condition` to your control condition (set to `"t0"` in this analysis) |

> 💡 **Tip:** Comments above each section explain what you might want to change. Look for text like "adjust according to your data" or "change according to your dataset".

### Output

Results are saved to `results/` in PDF and PNG formats:

```
results/
├── checkpoint_data/          # Cross-module checkpoints (see table above)
├── 01_quality_control/
├── 02_clustering_annotation/
├── 03_abundance/
├── 04_proximity/
└── 05_statistical_testing/
```

---

<br>

## Performed analysis description

Two analyses were run from this project. The **main analysis** (`modules/` → `results/`) covers all eight samples across both timepoints. A **product-only sub-analysis** (`modules_t0/` → `results_t0/`) repeats the workflow on the three `t0` CAR-T product samples in isolation. Both used a global seed of 42.

### Main analysis (all samples)

**Data loading and quality control.** Eight PXL files were loaded from `data/` according to `data/metadata.csv`, covering a 159-marker Proxiome Immuno panel (including the mIgG1, mIgG2a and mIgG2b isotype controls). Quality was assessed on read depth, average k-core, dangling-node fraction, UMI count per cell, isotype fraction and Tau score. Cells were retained if they were classified as `tau_type == "normal"`, had more than 20,000 UMIs (`n_umi_min_threshold = 20e3`; no upper bound, `n_umi_max_threshold = Inf`), and carried less than 0.1 % of their UMIs on isotype controls (`isotype_percent_threshold = 0.1`, applied as `isotype_fraction < 0.001`). Filtering retained **8,452 cells** (751–1,818 per sample).

**Marker selection.** Per-marker 50th and 95th CPM percentiles were computed and a 95th-percentile threshold of 2,000 CPM (`marker_cpm_p95_threshold = 2000`) was evaluated for diagnostic purposes, but the selection was overridden and **all 159 panel markers** were carried into dimensionality reduction (`selected_markers <- rownames(pg_data)`).

**Normalization and scaling.** Counts were normalized by centered log-ratio transformation across markers within each cell (`NormalizeData(normalization.method = "CLR", margin = 2)`) and scaled with `ScaleData()`.

**Dimensionality reduction and batch integration.** PCA was computed on all 159 markers (`RunPCA(npcs = 50)`), and UMAP on the first 10 principal components (`npcs = 10`, `n.neighbors = 30`). Samples were integrated with Harmony on `group.by.vars = "sample_alias"`, followed by UMAP on the first 20 Harmony components (`n.neighbors = 30`).

**Clustering and first-pass annotation.** A shared-nearest-neighbour graph was built on Harmony dimensions 1–10 and clustered with the Louvain algorithm at `resolution = 0.5`, yielding nine clusters. Clusters were manually annotated as CD4 T (clusters 0, 1), Mono (2), CD8 T (3, 6), Platelets (4, 7), NK (5) and DC (8).

**Annotation refinement by marker gating.** Because purely data-driven clustering recovered few T cells — particularly CD8⁺ and FMC63⁺ cells — in several samples, a guided strategy was adopted. Platelet clusters were removed, and the remaining cells were re-normalized (CLR, margin 2), re-scaled, re-embedded (PCA `npcs = 50`, UMAP on dimensions 1–10, `n.neighbors = 30`), re-integrated with Harmony on `group.by.vars = c("sample_alias", "condition")` with UMAP on Harmony dimensions 1–20, and re-clustered (Harmony dimensions 1–10, `resolution = 0.5`). NK, DC and Mono labels were carried over from the cluster annotation, with DC and Mono merged into a single **Myeloid** label. T cells were then gated on CLR values using `CD3e > 1.5`, `CD4 > 1`, `CD8 > 2.5` and `FMC63 > 2`, giving CD4 T FMC63−/+, CD8 T FMC63−/+, `T other` (CD3e⁺ but neither CD4⁺ nor CD8⁺) and `other`. The resulting gated labels (`manual_cell_gating`) were assigned to `cell_annotation` and drive modules 03–05. The final annotated object contains **7,159 cells**:

| Cell type | Cells |
|---|---:|
| CD4 T FMC63− | 2,065 |
| CD4 T FMC63+ | 760 |
| CD8 T FMC63− | 1,793 |
| CD8 T FMC63+ | 416 |
| Myeloid | 1,520 |
| NK | 440 |
| T other | 90 |
| other | 75 |

**Marker selection for statistical testing.** Background-associated markers were identified per cell type by partial least squares regression of marker abundance on the three isotype controls (`isotype_pls`, layer `"data"`, model matrix `~ factor(sample_alias)`, `remove_covariates = TRUE`), restricted to cell types present in at least two samples with at least two cells per sample. Markers were retained where mean CLR exceeded 10 isotype CLR standard deviations, or exceeded 3 standard deviations with a background PLS loading below 0.05. This gave 80–93 markers per cell type (`markers_to_test`).

**Protein abundance analysis.** Abundance distributions were generated for all 159 markers, on both CLR and raw count values, per sample, per cell type, per cell type faceted by sample, and per condition faceted by cell type (`abundance_plots_grouping_column = "sample_alias"`).

**Proximity (colocalization) analysis.** Proximity scores were summarized as the mean log2 ratio of observed over expected colocalization per marker pair, grouped by sample (`heatmap_grouping_column = "sample_alias"`). Heatmaps were restricted to the 60 most abundant markers of the `markers_to_test` union (`max_heatmap_markers = 60`), with marker ordering derived by hierarchical clustering of the first sample in the metadata (`donorD_t1`) and a colour scale capped at [−1, 1]. Equivalent heatmaps were produced per cell type and for a pooled T-cell annotation. Colocalization networks were drawn across all cell types and per cell type, with the FMC63-centred networks limited to the top 40 markers.

**Statistical testing.** For this module `condition` was redefined as `timepoint`, and `t0` was set as the baseline (`reference_condition = "t0"`), so all contrasts are `t1_vs_t0` (target − reference). Both statistical paths were enabled (`run_cell_level = TRUE`, `run_pseudobulk = TRUE`).

- *Cell-level.* Differential abundance was tested per cell type with `RunDAA` on that cell type's `markers_to_test`, using a Wilcoxon rank-sum test over pooled cells, followed by global FDR correction. Differential proximity was tested with `DifferentialProximityAnalysis` (`proximity_metric = "log2_ratio"`, `method = "seurat"`) for both self-proximity/clustering (`metric_type = "self"`) and colocalization (`metric_type = "all"`), with `min_cells_per_group = 10`, `diff_threshold = 0.01`, `min_pct = 0` and `min_diff_pct = -Inf`, and FDR correction applied after combining cell types. Differential abundance was additionally run restricted to **donorH** alone, using the same settings.
- *Pseudobulk.* CLR abundance was averaged per `sample_alias` × cell type, requiring at least 5 cells per sample within a cell type (`pseudobulk_min_cells = 5L`), and tested across samples with limma. The design was `~ 0 + condition` (one coefficient per condition, no intercept) with no additional covariates (`pseudobulk_covariate_cols = NULL`); contrasts were built numerically as target minus reference and fitted with `lmFit` → `contrasts.fit` → `eBayes` → `topTable`, followed by FDR correction. The same aggregation and design were applied to the clustering and colocalization proximity metrics.

Results were called significant at an adjusted p-value below 0.05 (`diff_p_adj_cutoff = 0.05`); cross–cell-type effect heatmaps additionally required an absolute effect size above 0.5 (`diff_abs_effect_min = 0.5`). Volcano plots label the 10 most significant markers per panel, and the pseudobulk summary colocalization heatmap shows the 50 strongest significant pairs.

### Product-only sub-analysis (`t0` samples)

The sample metadata was filtered to `timepoint == "t0"`, restricting the analysis to the three CAR-T product samples (donorF_t0, donorG_t0, donorH_t0).

**Quality control and marker selection.** The same thresholds were applied as in the main analysis (`tau_type == "normal"`, `n_umi > 20,000`, `n_umi_max_threshold = Inf`, isotype fraction below 0.1 %), retaining **2,369 cells** (807, 811 and 751 per sample). As above, `marker_cpm_p95_threshold = 2000` was evaluated but overridden, and all 159 markers were used.

**Normalization, dimensionality reduction and clustering.** Identical to the main analysis: CLR normalization (margin 2) and scaling; PCA with `npcs = 50` and UMAP on dimensions 1–10 with `n.neighbors = 30`; Harmony integration on `group.by.vars = "sample_alias"` with UMAP on Harmony dimensions 1–20; clustering on Harmony dimensions 1–10 at `resolution = 0.5`, with the same cluster-to-cell-type mapping (CD4 T: 0, 1; Mono: 2; CD8 T: 3, 6; Platelets: 4, 7; NK: 5; DC: 8).

**Annotation by marker gating.** No platelet-exclusion and reprocessing step was applied here. Gating was performed directly on all cells with cutoffs tuned to this subset — `CD3e > 1`, `CD4 > 1`, `CD8 > 2.75`, `FMC63 > 2` — and cells not matching any gate defaulted to `T other`. The annotated object contains **2,369 cells**: CD4 T FMC63− (1,000), CD4 T FMC63+ (687), CD8 T FMC63− (387), CD8 T FMC63+ (185) and T other (110). Isotype-PLS marker selection as described above gave 75–79 markers per cell type.

**Abundance and proximity.** Abundance distributions were generated for the first 10 panel markers plus FMC63, faceted by sample (`abundance_plots_grouping_column = "sample_alias"`). Colocalization heatmaps were grouped by sample with `max_heatmap_markers = 60` for the per-group heatmaps (59, 40 and 50 for the subsequent per-cluster, custom-marker and network heatmaps), using `donorF_t0` as the reference group for marker ordering.

**Statistical testing.** Module 05 was not run for this subset — the three `t0` samples share a single condition, `reference_condition` was left at the template placeholder, and no `results_t0/05_statistical_testing/` output exists.

### Software

Analyses were run under **R 4.4.3 (2025-02-28)** on macOS (aarch64), using pixelatorR 0.20.1, Seurat 5.4.0 (SeuratObject 5.3.0), harmony 2.0.4, limma 3.62.2, ComplexHeatmap 2.22.0, Matrix 1.7-5, ggraph 2.2.2, ggplotify 0.1.3, ggalluvial 0.12.5, future 1.70.0, tidyverse 2.0.0 (ggplot2 4.0.2, dplyr 1.2.0, tidyr 1.3.2, readr 2.1.5, purrr 1.2.1, tibble 3.3.1, stringr 1.6.0, forcats 1.0.0, lubridate 1.9.4) and here 1.0.2. Exact package versions are pinned in `renv.lock` at the root directory.

---

<br>

## Troubleshooting

<details>
<summary><strong>Package installation fails</strong></summary>

Run `renv::status()` and, if it reports differences, `renv::restore()` from the root directory.

If compilation fails, install compiler tools for your OS:
- **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS:** `xcode-select --install`
- **Linux:** `sudo apt install cmake libglpk-dev libhdf5-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev`

</details>

<details>
<summary><strong>"Cannot find function" error</strong></summary>

Run the Setup chunk in `modules/01_quality_control.qmd` first.

</details>

<details>
<summary><strong>File not found</strong></summary>

1. Open the project via `proxiome_analysis_template.Rproj`
2. Check that PXL files are in `data/`
3. Verify filenames in `metadata.csv` match exactly

</details>

<details>
<summary><strong>Checkpoint not found</strong></summary>

Run the upstream module first to create the missing file in `results/checkpoint_data/`:

- `pg_data_merged.rds`, `selected_markers.rds` → run module `01`
- `annotated_seurat_object.rds`, `markers_to_test.rds` → run module `02` to completion

</details>

<details>
<summary><strong>object 'pg_data' not found</strong></summary>

Run module `01` from the beginning in this R session, or run the first chunk of modules `03`–`05` to load from checkpoint.

</details>

<details>
<summary><strong>Memory errors</strong></summary>

- Close other applications
- Process fewer samples
- Use a machine with 16GB+ RAM

</details>

<br>

**Need more help?** [Open an issue](https://github.com/PixelgenTechnologies/proxiome-analysis-template/issues) · [pixelatorR issues](https://github.com/PixelgenTechnologies/pixelatorR/issues)

---

<br>

## Maintainers

Max Karlsson and Vincent van Hoef maintain the original Proxiome Analysis Template (PAT) code.

| | |
|---|---|
| **Max Karlsson** [@maxkarlsson](https://github.com/maxkarlsson) | [max.karlsson@pixelgen.com](mailto:max.karlsson@pixelgen.com) · [ORCID](https://orcid.org/0000-0002-7000-4416) |
| **Vincent van Hoef** [@vincent-van-hoef](https://github.com/vincent-van-hoef) | [vincent.vanhoef@pixelgen.com](mailto:vincent.vanhoef@pixelgen.com) · [ORCID](https://orcid.org/0000-0003-1707-7066) |

**Analysis performed by:** Lovisa Franzén (Field Bioinformatics Scientist, Pixelgen Technologies AB) — [lovisa.franzen@pixelgen.com](mailto:lovisa.franzen@pixelgen.com) · [@lfranzen](https://github.com/lfranzen) · [ORCID](https://orcid.org/0000-0003-3755-718X)

---

<br>

<div align="center">

**License:** [GNU GPL v2](LICENSE.md)

*Made with ❤️ by [Pixelgen Technologies](https://pixelgen.com)*

</div>
