# Set up conda environment and install packages

# Initialize env
conda env create -f environment.yml

# Activate env
conda activate suzuki_mpx_cart

# Update conda if needed
conda update -n base -c conda-forge conda

# Run within conda environment (in terminal)
conda install conda-forge::r-biocmanager

# conda install conda-forge::r-duckdb # Not working - install with `install.packages("duckdb")` in R instead
conda install conda-forge::r-conflicted conda-forge::r-readr conda-forge::r-purrr conda-forge::r-skimr

# conda install conda-forge::r-conflicted
# conda install conda-forge::r-readr
# conda install conda-forge::r-purrr
# conda install conda-forge::r-skimr

# conda install bioconda::r-harmony

# Visualizations pkgs
# conda install conda-forge::r-paletteer conda-forge::r-ggrepel conda-forge::r-ggplotify conda-forge::r-ggpubr bioconda::bioconductor-complexheatmap

conda install conda-forge::r-paletteer
conda install conda-forge::r-ggrepel
conda install conda-forge::r-ggplotify
conda install conda-forge::r-ggpubr

conda install bioconda::bioconductor-complexheatmap

# conda install conda-forge::r-gt

# VS Code specifics
# conda install conda-forge::r-httpgd r-languageserver
# conda install conda-forge::radian


# Install packages within R
R
install.packages("duckdb")
install.packages("harmony")
install.packages("arrow")
install.packages("gt")

remotes::install_github("PixelgenTechnologies/pixelatorR")

# BiocManager::install("ComplexHeatmap")
install.packages('tidyHeatmap', repos = c('https://stemangiola.r-universe.dev', 'https://cloud.r-project.org'))


BiocManager::install("pcaMethods") # Used for Seurat doublet prediction (PredictDoublets)

# Extra nice-to-have packages
install.packages('ggbeeswarm')
install.packages('ggplotify')
install.packages('ggpubr')


exit()

# Start R studio with conda env
/Applications/RStudio.app/Contents/MacOS/RStudio

# Start VS Code with conda env
/Applications/Visual\ Studio\ Code.app/Contents/MacOS/Electron

# Deactivate env
conda deactivate
