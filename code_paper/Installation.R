######################################################################
# This file is used to install the R package rome from 
# the author's Github repository: yk748/rome
######################################################################
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("yk748/rome")