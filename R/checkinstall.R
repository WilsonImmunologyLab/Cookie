if (!require("ggplot2")) install.packages("ggplot2")
if (!require("cluster")) install.packages("cluster")
if (!require("reshape2")) install.packages("reshape2")
if (!require("umap")) install.packages("umap")
if (!require("Rtsne")) install.packages("Rtsne")
if (!require("Rcpp")) install.packages("Rcpp")
if (!require("ggradar")) {
  if (!require("devtools")) {install.packages("devtools")}
  devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
}