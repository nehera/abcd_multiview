# This script is for installing packages onto MSI. 
# Note, the BIPnet package requires the external library gsl be available for installation.
# This can be done on MSI with the following command: module load gsl
r <- getOption("repos")
r["CRAN"] <- "http://cran.rstudio.com/"
options(repos=r)
install.packages("devtools")
library(devtools)
install_github('chekouo/BIPnet') 
install.packages("tidyverse")
install.packages("Rcpp")
install.packages("mvtnorm")