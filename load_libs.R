# TODO: Add comment
# 
# Author: slang
###############################################################################

## libpath <- '/home/slang/git_Projects/ExpressionSet/R/'
## source (paste(libpath,'../load_libs.R',sep=''))

#import(affy)  # Affymetrix pre-processing
#import(limma)
library(RSvgDevice)
library(gplots)
library(rgl)
#import(stats)
library(Rsubread)
library(DESeq) 
library("biomaRt")
library(stringr)
library(ggplot2)
library(reshape2)
library(scales)
library(boot)
source(paste(libpath,'ExpressionSet_Cluster.R', sep=''))
source(paste(libpath,'ExpressionSet.R', sep=''))
#source(paste(libpath,'ggplot_heatmap_fesability_study.R', sep=''))
source(paste(libpath,'NGSexpressionSet.R', sep=''))
source(paste(libpath,'PathwaysGMT.R', sep=''))
source(paste(libpath,'RFgrouping.R', sep=''))
source(paste(libpath,'SingleCells.R', sep=''))
source(paste(libpath,'Tool_RandomForest.R', sep=''))
