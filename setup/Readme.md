### Under construction

We will use R (version 4.0.5 or later) and Rstudio (version 1.4.1 or later) in this course.

Everybody should download Rstudio  (https://www.rstudio.com/) and install R (https://www.r-project.org/) before the course starts and the [required packages](Install_required_packages.R).

It is also preferable that you have all packages installed before the course starts (see below or [Install_required_packages.R](Install_required_packages.R). For those of you having trouble with installing the packages, we will try to find some time during the first computer lab to get everybody up and running.

Content of [required packages](Install_required_packages.R):

```R
# AB332_2021
# Please make sure that you have the required R packages installed.
# You can get a list of the package already installed on your computer by executing

installed.packages()

# in R alternatively you can just run the installation for each package
# to make sure that you have the latest version.
# Run these installation commands line-by-line in R (or Rstudio)
# and answer "yes" if you are asked to update any previously installed pakages:
# If you are asked if you want to compile from source it is usually easiest to answer "no"
# install.packages("ggplot2")   # To do plots, included in tidyverse package collection
install.packages("tidyverse") # To manipulate and visualize data https://www.tidyverse.org/packages/
install.packages("vegan")     # For ecology analyses
install.packages("ape")            # Phylogenetic tools
install.packages("recluster")      # Clustering tools
install.packages("dendextend")     # To work with dendrograms
install.packages("corrplot")       # makes nice correlation plots
install.packages("RcmdrMisc")      # diverse tools

# Packages from Bioconductor: (https://www.bioconductor.org/)
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(c("PCAtools","mixOmics"))

# Packages from github:
if (!requireNamespace("devtools", quietly = TRUE))  install.packages("devtools")
devtools::install_github("fawda123/ggord")
# This package might cause a problem. Here are some possible solutions:
# https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
# https://stackoverflow.com/questions/37776377/error-when-installing-an-r-package-from-github-could-not-find-build-tools-neces
```