library(pixelatorR)
library(dplyr)

##### COLORS #####

# Color palettes
cell_palette_unsorted <- c(
  PixelgenAccentColors(level = 5),
  PixelgenAccentColors(level = 8),
  PixelgenAccentColors(level = 3),
  PixelgenAccentColors(level = 10)
) |>
  as.character()

pal_samples <- c(
  "#D9E1ED",
  "#BECCE0",
  "#7B98C1",
  "#4B72A9",
  "#E3C2A0",
  "#C7A989",
  "#E8B6D5",
  "#CB73A8"
)

pal_time <- c(
  "#C0B2E9", #t0
  "#A1DAD0", #t1
  "#3FAC9B", #t2
  "#0D6D5E" #t3
)

heatmap_gradient <-
  colorRampPalette(
    c(
      "#1F395F",
      "#496389",
      "#728BB1",
      "#AABAD1",
      "#DFE5EE",
      "#FFFFFF",
      "#FFE0EA",
      "#E9AABF",
      "#CD6F8D",
      "#A23F5E",
      "#781534"
    )
  )

cherry_gradient <-
  colorRampPalette(
    c("gray90", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534")
  )

blues_gradient <-
  colorRampPalette(
    c("gray95", "#DFE5EE", "#AABAD1", "#728BB1", "#496389", "#1F395F")
  )

beige_gradient <-
  colorRampPalette(
    c("gray95", "#F0DCC6", "#D6B99A", "#C4A788", "#A68B6F", "#675440")
  )

# Bi-color palette for differential analysis
bi_col <- pixelatorR::Pixelgen_accent_colors |>
  dplyr::select(c("cyans", "pinks")) |>
  slice(7) |>
  unlist()

# Differential color palette
differential_pal <- colorRampPalette(c(bi_col[1], "white", bi_col[2]))

#' Annotated cell type color palette
cell_palette <-
  c(
    "CD4 T" = "#6D92D1",
    "CD4 Naive" = "#B9CDED",
    "Naive CD4 T" = "#B9CDED",
    "CD4 TSCM" = "#92B0E0",
    "CD4 TCM" = "#6D92D1",
    "CD4 TEM" = "#4A73C0",
    "CD4 TEFF" = "#224792",
    "CD4 TEMRA" = "#1C3A76",
    "Treg" = "#2955AE",
    "CD4 CTL" = "#2955AE",
    "CD8 T" = "#75C5B5",
    "CD8 Naive" = "#D0EDE6",
    "Naive CD8 T" = "#D0EDE6",
    "CD8 TSCM" = "#A2DACE",
    "CD8 TCM" = "#75C5B5",
    "CD8 TEM" = "#4AAF9D",
    "CD8 TEFF" = "#209785",
    "CD8 TEMRA" = "#1B7E6F",
    "CD4 T FMC63hi" = "#224792",
    # "CD4 T FMC63+" = "#224792",
    "CD4 T FMC63+" = "#171A6D",
    "CD4 T FMC63-" = "#5067D4",
    "CD8 T FMC63hi" = "#1B7E6F",
    # "CD8 T FMC63+" = "#1B7E6F",
    "CD8 T FMC63+" = "#40097A",
    "CD8 T FMC63-" = "#A75FC3",
    "other T" = "#3C7F96",
    "Other T" = "#3C7F96",
    "T other" = "#3C7F96",
    "MAIT" = "#0D6D5E",
    "DPT" = "#1B7E9F",
    "DNT" = "#7B787F",
    "dnT" = "#7B787F",
    # "NK" = "#A28EDB",
    "NK" = "#F5D56C",
    "CD56dim NK" = "#A28EDB",
    "CD56bright NK" = "#866CCD",
    "NK_CD56bright" = "#866CCD",
    "NKT" = "#6C4ABD",
    "B" = "#C7A989", # "#917557",
    "B naive" = "#F0DAC4",
    "Naive B" = "#F0DAC4",
    "B intermediate" = "#F0DAC4",
    "Memory B" = "#C7A989",
    "B memory" = "#C7A989",
    "Plasma cells" = "#DE9982",
    "Plasmablast" = "#DE9982",
    "DC" = "#DA94C1",
    "cDC1" = "#DA94C1",
    "cDC2" = "#BB5391",
    "pDC" = "#9C4579",
    # "Myeloid" = "#F7DFA0",
    "Myeloid" = "#ca1c87",
    "Mono" = "#F7DFA0",
    "Classical Mono" = "#F7DFA0",
    "CD14 Mono" = "#F7DFA0",
    "Intermediate Mono" = "#E6BE4A",
    "CD16 Mono" = "#BD9315",
    "Non-classical Mono" = "#BD9315",
    "other" = "#797979",
    "Neutrophils" = "#CDCDCD",
    "Basophils" = "#797979",
    "Platelet" = "#DB6365",
    "Platelets" = "#DB6365",
    "gdT" = "#5A3E9E",
    "discard" = "grey80"
  )


fmc_status_col <- c("grey80", "#CB73A8")
gradient_color <- c("#E5E5E5", "#A64463") # cherry_gradient(2)

rna_palette <- c(
  "CD4 T" = "#4968D9",
  "CD8 T" = "#AD5BCD",
  # "CD4 T FMC63hi" = "#224792",
  "CD4 T FMC63+" = "#161B6E",
  # "CD8 T FMC63hi" = "#1B7E6F",
  "CD8 T FMC63+" = "#5b3b74", # original color: "#4A0180"
  "Other T" = "#1B7E9F",
  "NK" = "#F9D849",
  "Myeloid" = "#CD7692",
  "Platelets" = "#797979"
)

cell_type_palette_fmc <-
  c(
    "CD4 T" = "#5067D4",
    "CD4 T FMC63-" = "#5067D4",
    "CD4 T FMC63+" = "#171A6D",
    "CD8 T" = "#A75FC3",
    "CD8 T FMC63-" = "#A75FC3",
    "CD8 T FMC63+" = "#5b3b74", # original color: "#4A0180"
    "DC" = "#DA94C1",
    "Mono" = "#C47892",
    "Myeloid" = "#ca1c87",
    "NK" = "#F5D56C",
    "Platelets" = "#DB6365",
    "Other T" = "#3C7F96",
    "T other" = "#3C7F96"
  )


pal_org_bl <- c("#053061", "white", "#670020")

colors <- c(
  "#053061",
  "#2166AC",
  "#4393C3",
  "#92C5DE",
  "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7",
  "#F4A582",
  "#D6604D",
  "#B2182B",
  "#67001F"
)

cherry_gradient <-
  colorRampPalette(
    c("gray90", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534")
  )

retro_palette <- colorRampPalette(c("#77BDDD", "#fcf4ad", "#FB2340"))

pal_patient <- c(
  "Patient 1" = "#7B98C1",
  "Patient 2" = "#C7A989",
  "Patient 3" = "#CB73A8"
)

# pal_donor <- c( #old
#   donorA = "#7B98C1",
#   donorB = "#C7A989",
#   donorC = "#CB73A8",
#   donorD = "#466EC0",
#   donorE = "#3FAC9B",
#   donorF = "#cd6b41",
#   donorG = "#6C4ABD",
#   donorH = "#773031"
# )

pal_donor <- c(
  # responders: blues, teals and violets, separated mainly by lightness
  donorA = "#1A4568", # dark steel blue
  donorB = "#007476", # dark teal
  donorD = "#74BDBF", # light teal
  donorE = "#AFF4F7", # pale cyan
  donorF = "#3B0AAC", # deep indigo
  donorG = "#704D9E", # muted purple
  donorH = "#A25BD7", # medium violet
  donorJ = "#A381D5", # soft violet
  donorK = "#D5AAFE", # light lavender
  donorL = "#7287AC", # slate blue
  donorM = "#4A007E", # dark purple
  # non-responders: brown, orange, yellow
  donorC = "#954A0E", # brown
  donorI = "#E8710A", # orange
  donorN = "#F5C242" # yellow
)

##### THEMES #####
# Custom ggplot themes

theme_violin <- function(base_size = 12) {
  theme_linedraw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid = element_blank()
    )
}

theme_umap <- function(legend.pos = "right") {
  theme(
    legend.position = legend.pos,
    aspect.ratio = 1,
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)
  )
}

##### HELPER FUNCTIONS #####
library(DT)
display_dt <- function(x, ...) {
  DT::datatable(
    x,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
      pageLength = 15
    ),
    ...
  )
}
