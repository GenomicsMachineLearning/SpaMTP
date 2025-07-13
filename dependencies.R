#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://mirror.aarnet.edu.au/pub/CRAN/"))
# Basic packages - remote installation and load/save packages
install.packages("remotes", dependencies = FALSE)

install.packages("tiff", version = "0.1.12", upgrade = FALSE, dependencies = FALSE)
install.packages("fftwtools", version = "0.9.11", upgrade = FALSE, dependencies = FALSE)

# Version 4.47.1
remotes::install_github("aoles/EBImage", upgrade = FALSE, dependencies = FALSE)

# Latest/Version 2.11.1
remotes::install_github("kuwisdelu/matter", upgrade = FALSE, dependencies = FALSE)
# Latest/Version 1.5.0
remotes::install_github("kuwisdelu/CardinalIO", upgrade = FALSE, dependencies = FALSE)
# Latest/Version 3.9.0
remotes::install_github("kuwisdelu/Cardinal", upgrade = FALSE, dependencies = FALSE)

# Latest/Version 1.7.0
remotes::install_github("lifs-tools/rgoslin", upgrade = FALSE, dependencies = FALSE)
# Version 3.21
remotes::install_github("MarioniLab/DropletUtils@RELEASE_3_21", upgrade = FALSE, dependencies = FALSE)

# Latest
remotes::install_github("GenomicsMachineLearning/SpaMTP", upgrade = FALSE, dependencies = FALSE)
