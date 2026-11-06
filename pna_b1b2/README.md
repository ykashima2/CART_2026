<div align="center">

# Proxiome Analysis Template

**A complete, reproducible workflow for single-cell Proximity Network Analysis**

[![R](https://img.shields.io/badge/R-4.4.3-blue?logo=r)](https://www.r-project.org/)
[![Quarto](https://img.shields.io/badge/Quarto-Document-purple)](https://quarto.org/)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-green.svg)](../pna_1/LICENSE.md)

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
- [Troubleshooting](#troubleshooting)
- [Maintainers](#maintainers)

---

<br>

## Sample Information

This project analyses **both PNA batches together**: **16 samples** from **11 donors**, covering **two CAR constructs**. Samples were collected either from the CAR-T product itself (`t0`) or from PBMCs drawn after infusion (`t1`).

The PXL files are **read in place** from the two source projects — this project's `data/` folder holds only `metadata.csv`.

| sample_id | sample_alias | donor | timepoint | dpi | sample_type | cell_viability (%) | construct | status | batch | condition |
|-----------|--------------|-------|:---------:|:---:|-----------|:------------------:|-----------|--------|:-----:|----------|
| TKB_CAR-T_10 | donorD_t1 | donorD | t1 | 7 | PBMC | 76 | tisa-cel | responding | b1 | t1_tisa-cel |
| TKB_CAR-T_2 | donorE_t1 | donorE | t1 | 7 | PBMC | 59 | tisa-cel | responding | b1 | t1_tisa-cel |
| TKB_CAR-T_50 | donorF_t0 | donorF | t0 | 0 | CAR-T_PRODUCT | 9 | tisa-cel | responding | b1 | t0_tisa-cel |
| TKB_CAR-T_55 | donorF_t1 | donorF | t1 | 11 | PBMC | 95 | tisa-cel | responding | b1 | t1_tisa-cel |
| TKB_CAR-T_194 | donorG_t0 | donorG | t0 | 0 | CAR-T_PRODUCT | 29 | tisa-cel | responding | b1 | t0_tisa-cel |
| TKB_CAR-T_198 | donorG_t1 | donorG | t1 | 7 | PBMC | 87 | tisa-cel | responding | b1 | t1_tisa-cel |
| TKB_CAR-T_175 | donorH_t0 | donorH | t0 | 0 | CAR-T_PRODUCT | 37 | tisa-cel | responding | b1 | t0_tisa-cel |
| TKB_CAR-T_181 | donorH_t1 | donorH | t1 | 10 | PBMC | 96 | tisa-cel | responding | b1 | t1_tisa-cel |
| TKB_CAR-T_217 | donorI_t1 | donorI | t1 | NA | PBMC | NA | axi-cel | non-responding | b2 | t1_axi-cel |
| TKB_CAR-T_227 | donorJ_t1 | donorJ | t1 | NA | PBMC | NA | axi-cel | responding | b2 | t1_axi-cel |
| TKB_CAR-T_177 | donorK_t0 | donorK | t0 | NA | CAR-T_PRODUCT | NA | axi-cel | responding | b2 | t0_axi-cel |
| TKB_CAR-T_180 | donorK_t1 | donorK | t1 | NA | PBMC | NA | axi-cel | responding | b2 | t1_axi-cel |
| TKB_CAR-T_212 | donorL_t0 | donorL | t0 | NA | CAR-T_PRODUCT | NA | axi-cel | responding | b2 | t0_axi-cel |
| TKB_CAR-T_221 | donorL_t1 | donorL | t1 | NA | PBMC | NA | axi-cel | responding | b2 | t1_axi-cel |
| TKB_CAR-T_278 | donorM_t1 | donorM | t1 | NA | PBMC | NA | axi-cel | responding | b2 | t1_axi-cel |
| TKB_CAR-T_301 | donorN_t1 | donorN | t1 | NA | PBMC | NA | axi-cel | non-responding | b2 | t1_axi-cel |

**At a glance**

| | |
|---|---|
| Samples | 16 |
| Donors | 11 (donorD–donorN) |
| Batches | `b1` (8, from `pna_1`) · `b2` (8, from `pna_2`) |
| Constructs | `tisa-cel` (8, = b1) · `axi-cel` (8, = b2) |
| Sample types | `CAR-T_PRODUCT` (5) · `PBMC` (11) |
| Timepoints | `t0` (5) · `t1` (11) |
| Conditions | `t0_tisa-cel` (3) · `t1_tisa-cel` (5) · `t0_axi-cel` (2) · `t1_axi-cel` (6) |
| Paired donors (t0 + t1) | donorF, donorG, donorH, donorK, donorL |
| Status | responding (14) · non-responding (2: donorI, donorN) |

> 📁 **Where the data lives.** `data/` contains only `metadata.csv`. Each `file_path` is relative to the root directory and points at `pna_1/data/<sample_id>.layout.pxl` or `pna_2/data/<sample_id>.layout.pxl`. Module `01` resolves them with `here::here("..", metadata$file_path)`. Do not move, rename or delete `pna_1/` or `pna_2/`. The `.pxl` files are not included in this repository: to replicate this analysis, place the PXL files listed in `data/metadata.csv` in `pna_1/data/` and `pna_2/data/`.

> ⚠️ Do not use `dpi` or `cell_viability` as limma covariates — they are not available for all samples, so `model.matrix` would drop those samples.

> ⚠️ **Batch is confounded with construct.** `batch` and `construct` are the same split (b1 ≡ tisa-cel, b2 ≡ axi-cel), and every donor belongs to exactly one construct. Any batch effect is therefore inseparable from a construct effect, and adding `batch` to `pseudobulk_covariate_cols` in module `05` would make the design collinear. See that module's header comments.

> 📝 This table mirrors `data/metadata.csv`. If you edit the CSV (samples, order, values), update this section to match. Row order sets factor levels and therefore plot order.

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

This project is already set up. Unlike a standalone template copy, it holds **no PXL files of its own** — `data/` contains only `metadata.csv`, and each row points at a `.pxl` file in the sibling `pna_1/` or `pna_2/` project.

### Project structure

```
root directory/
├── pna_1/data/*.layout.pxl                ← Batch 1 PXL files (read in place, 8 samples)
├── pna_2/data/*.layout.pxl                ← Batch 2 PXL files (read in place, 8 samples)
└── pna_b1b2/                              ← THIS project
    ├── data/
    │   └── metadata.csv                   ← Sample info for all 16 samples (no PXL files here)
    ├── modules/                           ← Analysis notebooks (run in order)
    │   ├── common_setup.R                 ← Shared setup (sourced automatically)
    │   ├── marker_selection_functions.R   ← CLR-based marker selection helpers
    │   ├── platelet_cleanup_functions.R   ← Helpers for module 02b
    │   ├── 01_quality_control.qmd         ← Start here (load data + QC)
    │   ├── 02_clustering_annotation.qmd
    │   ├── 02b_platelet_cleanup.qmd       ← Platelet removal from T cell graphs (long run)
    │   ├── 03_abundance.qmd
    │   ├── 04_proximity.qmd
    │   ├── 05_statistical_testing.qmd     ← Optional
    │   ├── 05b_statistical_testing_response.qmd ← Optional (responding vs non-responding, axi-cel t1)
    │   ├── 06_cell_visualization.qmd      ← Optional
    │   └── 08_manuscript_figures.qmd      ← Figure 7 and supplementary panels
    ├── modules_cleaned/                   ← Modules 03, 04, 05, 05b on the platelet-cleaned T cell object
    │   ├── 00_setup_cleaned.R
    │   ├── 03_abundance_cleaned.qmd
    │   ├── 04_proximity_cleaned.qmd
    │   ├── 05_statistical_testing_cleaned.qmd
    │   └── 05b_statistical_testing_response_cleaned.qmd
    ├── results/                           ← Created automatically when running the modules
    ├── results_cleaned/                   ← Output of modules_cleaned/
    ├── proxiome_analysis_template.Rproj   ← Double-click to open the project in RStudio
    ├── proxiome_analysis_template.qmd     ← Full re-render of modules 01–05b (optional)
    └── proxiome_analysis_platelet_cleaned.qmd ← Full re-render of the cleaned modules (optional)
```

### Create metadata.csv

The `metadata.csv` file describes the samples and links them to the corresponding PXL files. 

| Column | Required | Description | Example |
|--------|:--------:|-------------|---------|
| `sample_id` | ✓ | Unique identifier | `S1` |
| `sample_alias` | ✓ | Descriptive name | `donorD_t1` |
| `condition` | ✓ | `timepoint`_`construct`; the grouping variable for the contrasts in module `05` | `t0_axi-cel` |
| `file_path` | ✓ | Path to the PXL file, **relative to the root directory** | `pna_1/data/TKB_CAR-T_10.layout.pxl` |
| `batch` | | Source project: `b1` = `pna_1`, `b2` = `pna_2` | `b1` |
| `sample_type` | | Sample origin: `CAR-T_PRODUCT` (t0) or `PBMC` (t1) | `PBMC` |
| *other* | | Any extra columns (`donor`, `timepoint`, etc.) | - |

> 💡 Extra columns are available as variables in the template for grouping or covariates.

**Example rows** (the first of each batch — see [Sample Information](#sample-information) for the full 16-row table):

```csv
sample_id,sample_alt_name,sample_alias,dpi,donor,cell_viability,sample_type,file_path,timepoint,construct,status,batch,condition
TKB_CAR-T_10,TKB_CAR-T_10,donorD_t1,7,donorD,76,PBMC,pna_1/data/TKB_CAR-T_10.layout.pxl,t1,tisa-cel,responding,b1,t1_tisa-cel
TKB_CAR-T_217,TKB_CAR-T_217,donorI_t1,NA,donorI,NA,PBMC,pna_2/data/TKB_CAR-T_217.layout.pxl,t1,axi-cel,non-responding,b2,t1_axi-cel
```

> 📝 Row order controls plot order. Both `,` and `;` separators work. Every `file_path` is resolved as `here::here("..", file_path)` from this project, so it must exist under the root directory.

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
4. **Optional:** Open `modules/05_statistical_testing.qmd` for differential testing, `modules/05b_statistical_testing_response.qmd` for responding vs non-responding testing (axi-cel, t1) and `modules/06_cell_visualization.qmd` for single-cell graph plots

> 💡 **Tip:** Each module starts with a short intro. Search for **Decision checkpoint** to find values you need to edit.

### Iterate on downstream modules

After module `02` finishes, you can open `03`, `04`, `05`, or `06` in a **new R session** and run the first chunk — packages, palettes, and saved data load automatically.

| Checkpoint file | Created by | Used by |
|-----------------|------------|---------|
| `pg_data_merged.rds` | module `01` (filtered cells) | module `02` |
| `selected_markers.rds` | module `01` | module `02` |
| `annotated_seurat_object.rds` | module `02` | modules `03`–`06` |
| `markers_to_test.rds` | module `02` | modules `04`–`05` |
| `platelet_tcell_protein_weights.rds` | module `02` | module `02b` |
| `pg_data_platelet_cleaned.rds` | module `02b` | `modules_cleaned/`, module `08` |

All checkpoints are saved under `results/checkpoint_data/`.

### Platelet-cleaned workstream

Platelet contamination in the T cells is removed at the cell-graph level in module `02b`, and the T cell analyses are then repeated on the cleaned object:

1. Run modules `01`–`02` (module `02` saves the platelet / T cell labels and protein weights).
2. Run `modules/02b_platelet_cleanup.qmd` with `run_cleanup <- TRUE`. This is a long run with high memory use; it is done in resumable batches and writes `pg_data_platelet_cleaned.rds` (T cells only, recomputed counts and proximity scores).
3. Run the modules in `modules_cleaned/` (`03`, `04`, `05`, `05b`; or render `proxiome_analysis_platelet_cleaned.qmd`). They read the cleaned object and write to `results_cleaned/`.
4. `modules/08_manuscript_figures.qmd` generates the Figure 7 and supplementary panels.

`modules_cleaned/04_proximity_cleaned.qmd` also contains the CD3e-centred colocalization networks (per construct and timepoint) and the FMC63 clustering-score heatmap at t0.

### Related project

`../pna_3` is a separate axi-cel validation set analysed with the same workflow. Its module `08` (figS15–S20) reads this project's platelet-cleaned object (`results/checkpoint_data/pg_data_platelet_cleaned.rds`) read-only, so keep this folder's name and location unchanged.

### Full re-render (optional)

To regenerate one HTML report of the entire workflow, install [Quarto](https://quarto.org/) and render `proxiome_analysis_template.qmd`. For day-to-day analysis, use the individual modules instead.

### ⚠️ Where you need to make choices

| Module | What to do |
|--------|------------|
| **01 — sample table** | Check `metadata.csv` paths and sample IDs match your data |
| **01 — QC thresholds** | Set `n_umi_min_threshold`, `n_umi_max_threshold`, `isotype_percent_threshold` after reviewing QC plots; re-run section 1.3.5 onward |
| **01 — marker selection** | Review `selected_markers` after running the selection code |
| **02 — manual annotation** | **Critical step!** Edit `cluster_cell_annotation` to assign cell type names to clusters |
| **05 — reference condition** (optional) | Contrasts are set as explicit (target, reference) pairs in `contrasts_to_test`; `reference_condition <- "t0_tisa-cel"` sets the limma baseline; `pseudobulk_covariate_cols <- c("donor")` |
| **05b — response subset** (optional) | Tests `status` (responding vs non-responding, `reference_condition <- "non-responding"`) on axi-cel t1 samples only; set `response_timepoint`, `response_construct` and `response_column` to change the subset. 4 vs 2 donors, so results are exploratory |
| **02b — clean-up settings** | Set `max_contam_fraction`, `memory_budget_gb` and the segmentation settings; set `run_cleanup <- TRUE` to start the batch run (long) |
| **04 cleaned — CD3e network / clustering heatmap** (optional) | Marker and cell-type settings are at the top of each chunk (`marker_focus`, `partner_markers`, `net_celltype`, `clust_celltype`, abundance and join-count thresholds) |
| **08 — manuscript figures** (optional) | Generates the Figure 7 and supplementary panels; needs the platelet-cleaned object (module `02b`) and the MPX results in `mpx_2/results/b1b2_t0t1/` |
| **06 — cell type & markers** (optional) | Set `cell_type_of_interest` and `markers_to_viz` to a cell type and markers that cluster/colocalize in your data |

> 💡 **Tip:** Comments above each section explain what you might want to change. Look for text like "adjust according to your data" or "change according to your dataset".

### Output

Results are saved to `results/` in PDF and PNG formats:

```
results/
├── checkpoint_data/          # Cross-module checkpoints (see table above)
├── 01_quality_control/
├── 02_clustering_annotation/
├── 02b_platelet_cleanup/     # Batch files (counts, proximity scores, QC, cell graphs)
├── 03_abundance/
├── 04_proximity/
├── 05_statistical_testing/
├── 05b_statistical_testing_response/
├── 06_cell_visualization/
└── 08_manuscript_figures/    # figure7/ and supplement/ panels (png + pdf)

results_cleaned/              # Output of modules_cleaned/
├── 03_abundance/
├── 04_proximity/
├── 05_statistical_testing/
└── 05b_statistical_testing_response/
```

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
2. Check that the PXL files are still in `../pna_1/data/` and `../pna_2/data/` — this project reads them in place
3. Verify that each `file_path` in `metadata.csv` resolves from the root directory (e.g. `pna_1/data/TKB_CAR-T_10.layout.pxl`)

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

**License:** [GNU GPL v2](../pna_1/LICENSE.md)

*Made with ❤️ by [Pixelgen Technologies](https://pixelgen.com)*

</div>
