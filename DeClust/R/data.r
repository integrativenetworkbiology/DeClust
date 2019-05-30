
#' Immune and stromal markers from R package estimate.
#'
#' @format A data frame with 2 rows and 142 columns:
"SI_geneset"



#' bulky gene expression for testing
#' @format  A matrix with 200 rows(genes) and 102 columns(samples)
"exprM"

#' precomputed subtypes for 13 TCGA datasets by DeClust
#' @format  A list of 13 items. Each item corresponds to one TCGA dataset, and is a vector of subtype annotation for samples in that TCGA dataset. 
"TCGAdeClustsubtype"

#' precomputed expression profiles for 13 TCGA datasets by DeClust
#'
#'DeClust was run for each of the 13 TCGA datasets separately (see methods of the paper for details). After that, the estimated profiles for cancer subtypes and immune/stromal compartments were assembled together into one single matrix for down-stream analysis. Only genes existed in more than half of the 13 datasets were kept. In order to make profiles more comparable between different datasets (e.g.,removing batch effects), quantile normalization was applied to the assemebled profiles.

#' @format  A matrix of 10890 rows(genes) and 87 columns(different cancer subtypes,and tissue-specific immune and stromal compartments) 
"TCGAdeClustprofileM"



