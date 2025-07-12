#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://mirror.aarnet.edu.au/pub/CRAN/"))
# Basic packages - remote installation and load/save packages
install.packages("remotes", dependencies = FALSE)

install.packages("tiff", upgrade = FALSE, dependencies = FALSE)
install.packages("fftwtools", upgrade = FALSE, dependencies = FALSE)

remotes::install_github("aoles/EBImage", upgrade = FALSE, dependencies = FALSE)

remotes::install_github("kuwisdelu/matter", upgrade = FALSE, dependencies = FALSE)
remotes::install_github("kuwisdelu/CardinalIO", upgrade = FALSE, dependencies = FALSE)
remotes::install_github("kuwisdelu/Cardinal", upgrade = FALSE, dependencies = FALSE)

remotes::install_github("lifs-tools/rgoslin", upgrade = FALSE, dependencies = FALSE)
remotes::install_github("MarioniLab/DropletUtils@RELEASE_3_21", upgrade = FALSE, dependencies = FALSE)

remotes::install_github("GenomicsMachineLearning/SpaMTP", upgrade = FALSE, dependencies = FALSE)
