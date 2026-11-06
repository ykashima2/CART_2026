# pna_3: axi-cel validation set

Stand-alone validation analysis of five axi-cel PNA samples, using the same workflow as `../pna_b1b2` (QC, Harmony clustering, NMF reference annotation, platelet clean-up, abundance, proximity). Results are meant to be compared with the axi-cel samples of `pna_b1b2`.

| sample_alias | sample_id | timepoint | sample_type |
|:--|:--|:--|:--|
| donorO_t0 | TKB_CAR-T_327 | t0 | CAR-T_PRODUCT |
| donorO_t1 | TKB_CAR-T_328 | t1 | PBMC |
| donorP_t1 | TKB_CAR-T_281 | t1 | PBMC |
| donorQ_t1 | TKB_CAR-T_300 | t1 | PBMC |
| donorR_t0 | TKB_CAR-T_310 | t0 | CAR-T_PRODUCT |

Only donorO is paired, so the t1 vs t0 test in module 05 (the only contrast) is exploratory.

> 📁 **PXL files.** The `.pxl` files are not included in this repository. To replicate this analysis, place the PXL files listed in `data/metadata.csv` (`pna_3/data/<sample_id>.layout.pxl`) in the `data/` folder of this project.

R packages come from the shared renv project at the root directory (`renv::restore()`; versions in `renv.lock`). Open the project through `proxiome_analysis_template.Rproj`.

## Run order

`modules/01` -> `02` -> `03` -> `04` -> `05` (all via `proxiome_analysis_template.qmd`), then optionally `02b` (platelet clean-up, set `run_cleanup <- TRUE`; long run, high memory), then `modules_cleaned/` (`proxiome_analysis_platelet_cleaned.qmd`) and `08_manuscript_figures.qmd`.

Differences from `pna_b1b2`:

- `data/metadata.csv`: `condition` is `timepoint_construct` (as in b1b2); the sample type is in `sample_type`.
- Module 02: platelet clusters are picked from the data (clusters where >= 50% of cells are labelled Platelets), not fixed.
- Module 05: single contrast `t1_axi-cel` vs `t0_axi-cel`, no donor blocking.
- Module 08: supplementary figures `figS1`-`figS14` only, no differential results; S10-S14 show raw t0 vs t1 scores. Figures are written to `results/08_manuscript_figures/`. Uses the platelet-cleaned object when present, otherwise `pg_data`.
- Module 08, section 8.12: `figS15`-`figS20` compare pna_3 with the axi-cel samples of `pna_b1b2` (read-only; FMC63 colocalization heatmaps per batch and timepoint and per sample, concordance scatter, selected pairs, raw abundance, sample PCA). They only run once the platelet-cleaned objects of both projects exist, and are skipped otherwise.
- No 05b, 06 or 07.

## Maintainers

The workflow is based on the Proxiome Analysis Template (PAT), maintained by Max Karlsson and Vincent van Hoef.

**Analysis performed by:** Lovisa Franzén (Field Bioinformatics Scientist, Pixelgen Technologies AB) — [lovisa.franzen@pixelgen.com](mailto:lovisa.franzen@pixelgen.com) · [@lfranzen](https://github.com/lfranzen) · [ORCID](https://orcid.org/0000-0003-3755-718X)
