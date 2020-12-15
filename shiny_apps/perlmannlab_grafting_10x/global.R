## Written 2017 by Asa Bjorklund asa.bjorklund@scilifelab.se

# All files created with script analysis/grafting/bin/extract_data_for_shiny_v2.R

library(shiny)
library(data.table)
library(ggplot2)
library(gridExtra)
library(Matrix)

load("integrated_rat_scale_res0.1.Rdata")
# contains: meta,expr,genes,


