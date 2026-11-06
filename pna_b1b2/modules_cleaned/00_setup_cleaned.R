# Setup for the platelet-cleaned workstream (modules 03, 04, 05, 05b).
# Sources the original common_setup.R unchanged, then points the cleaned copies
# at the cleaned object and the separate output folder.

source(here::here("pna_b1b2", "modules", "common_setup.R"))

# Output root for this workstream (module copies write to results_cleaned/NN_*)
results_root <- here::here("pna_b1b2", "results_cleaned")

# Cleaned object written by module 02b (T cells only, platelet-cleaned counts,
# recomputed proximity scores held in memory rather than in a PXL database)
input_object_name <- "pg_data_platelet_cleaned.rds"

# Cell types present in the cleaned object
tcell_types <- c("CD4 T FMC63+", "CD4 T FMC63-", "CD8 T FMC63+", "CD8 T FMC63-")

# Marker lists per cell type, restricted to the T-cell types above
load_markers_to_test_cleaned <- function() {
  load_checkpoint("markers_to_test.rds")[tcell_types]
}

# Guard against a markers_to_test already in the session from the original run
if (exists("markers_to_test", inherits = TRUE)) {
  markers_to_test <- markers_to_test[intersect(names(markers_to_test), tcell_types)]
}
