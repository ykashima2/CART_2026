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
- [Performed analysis description](#performed-analysis-description)
- [Troubleshooting](#troubleshooting)
- [Maintainers](#maintainers)

---

<br>

## Sample Information

This analysis covers **8 PNA samples** from **6 donors**, all carrying the **axi-cel** CAR construct. Samples were collected either from the CAR-T product itself (`t0`) or from PBMCs drawn after infusion (`t1`).

| sample_id | sample_alias | donor | timepoint | dpi | condition | cell_viability (%) | construct | status |
|-----------|--------------|-------|:---------:|:---:|-----------|:------------------:|-----------|--------|
| TKB_CAR-T_217 | donorI_t1 | donorI | t1 | NA | PBMC | NA | axi-cel | non-responding |
| TKB_CAR-T_227 | donorJ_t1 | donorJ | t1 | NA | PBMC | NA | axi-cel | responding |
| TKB_CAR-T_177 | donorK_t0 | donorK | t0 | NA | CAR-T_PRODUCT | NA | axi-cel | responding |
| TKB_CAR-T_180 | donorK_t1 | donorK | t1 | NA | PBMC | NA | axi-cel | responding |
| TKB_CAR-T_212 | donorL_t0 | donorL | t0 | NA | CAR-T_PRODUCT | NA | axi-cel | responding |
| TKB_CAR-T_221 | donorL_t1 | donorL | t1 | NA | PBMC | NA | axi-cel | responding |
| TKB_CAR-T_278 | donorM_t1 | donorM | t1 | NA | PBMC | NA | axi-cel | responding |
| TKB_CAR-T_301 | donorN_t1 | donorN | t1 | NA | PBMC | NA | axi-cel | non-responding |

**At a glance**

| | |
|---|---|
| Samples | 8 |
| Donors | 6 (donorI–donorN) |
| Conditions | `CAR-T_PRODUCT` (2) · `PBMC` (6) |
| Timepoints | `t0` (2) · `t1` (6) |
| Paired donors (t0 + t1) | donorK, donorL |
| Construct | axi-cel |
| Status | responding (6: donorJ–M) · non-responding (2: donorI, donorN) |

> 📁 **PXL files.** The `.pxl` files are not included in this repository. To replicate this analysis, place the PXL files listed in `data/metadata.csv` (`<sample_id>.layout.pxl`) in the `data/` folder of this project.

> 📝 This table mirrors `data/metadata.csv`. If you edit the CSV (samples, order, values), update this section to match. Each `sample_id` corresponds to `<sample_id>.layout.pxl` in `data/`.

> 💡 This is the **axi-cel** cohort. The parallel **tisa-cel** cohort lives in `pna_1/`; both metadata tables share an identical 11-column schema so the two can be compared by `construct`.

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
pna_2/
├── data/
│   ├── metadata.csv                  ← Your sample info
│   └── *.layout.pxl                  ← Your PXL files
├── modules/                          ← Analysis notebooks (run in order)
│   ├── common_setup.R                ← Shared setup (sourced automatically)
│   ├── 01_quality_control.qmd        ← Start here (load data + QC)
│   ├── 02_clustering_annotation.qmd
│   ├── 03_abundance.qmd
│   ├── 04_proximity.qmd
│   └── 05_statistical_testing.qmd    ← Optional
├── results/                          ← Created automatically when running the PAT
├── proxiome_analysis_template.Rproj  ← Double-click to open the project in RStudio
└── proxiome_analysis_template.qmd  ← Full re-render of all modules (optional)
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

Add one row per sample. The first two rows of `data/metadata.csv` in this project (see [Sample Information](#sample-information) for the full table):

```csv
sample_id,sample_alt_name,sample_alias,dpi,donor,cell_viability,condition,file_path,timepoint,construct,status
TKB_CAR-T_217,TKB_CAR-T_217,donorI_t1,NA,donorI,NA,PBMC,TKB_CAR-T_217.layout.pxl,t1,axi-cel,non-responding
TKB_CAR-T_227,TKB_CAR-T_227,donorJ_t1,NA,donorJ,NA,PBMC,TKB_CAR-T_227.layout.pxl,t1,axi-cel,responding
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

The analysis covers the eight axi-cel PXL files in `data/` (donors I–N), loaded according to `data/metadata.csv`. A global seed of 42 was used.

**Quality control.** Cells were retained if they were classified as `tau_type == "normal"`, had more than 20,000 UMIs (`n_umi_min_threshold = 20e3`; no upper bound) and carried less than 0.1 % of their UMIs on isotype controls (`isotype_percent_threshold = 0.1`).

**Marker selection.** A 95th-percentile threshold of 2,000 CPM (`marker_cpm_p95_threshold = 2000`) was evaluated, but all panel markers were carried into dimensionality reduction (`selected_markers <- rownames(pg_data)`).

**Normalization, integration and clustering.** Counts were CLR-normalized across markers within each cell and scaled. PCA used 50 components and UMAP the first 10 (`n.neighbors = 30`). Samples were integrated with Harmony on `sample_alias`, followed by UMAP on the first 20 Harmony components. Louvain clustering on Harmony dimensions 1–10 at `resolution = 0.5` gave nine clusters, annotated as Platelets (0, 5), CD4 T (1), Myeloid (2, 7, 8), CD8 T (3, 4) and T other (6).

**Annotation refinement by marker gating.** Platelet clusters were removed and the remaining cells were re-embedded (same settings) and re-integrated with Harmony on `sample_alias` and `condition`. T cells were gated on CLR values using `CD3e > 1.5`, `CD4 > 1` and `CD8 > 2`. FMC63+ cells were called with `FMC63 > 1.75` at t0 and `FMC63 > 0.5` at t1.

**Marker selection for statistical testing.** Background-associated markers were identified per cell type by partial least squares regression on the isotype controls (`isotype_pls`). Markers were retained where mean CLR exceeded 10 isotype CLR standard deviations, or exceeded 3 standard deviations with a background PLS loading below 0.05 (`markers_to_test`).

**Abundance and proximity.** Abundance distributions and colocalization heatmaps were grouped by `sample_alias`. Heatmaps showed the 60 most abundant markers, and 40 in the later per-cluster heatmaps.

**Statistical testing.** `t0` was the reference condition (contrast `t1` vs `t0`). Cell-level tests and limma pseudobulk (`~ 0 + condition`, no covariates, at least 5 cells per sample and cell type) were both run. Proximity tests used `log2_ratio`, `method = "seurat"`, `min_cells_per_group = 10` and `diff_threshold = 0.01`. Results were called significant at adjusted p < 0.05, and cross-cell-type heatmaps also required an absolute effect above 0.5.

Exact package versions are pinned in `renv.lock` at the root directory.

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

**License:** [GNU GPL v2](../pna_1/LICENSE.md)

*Made with ❤️ by [Pixelgen Technologies](https://pixelgen.com)*

</div>
