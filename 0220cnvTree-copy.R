library(cnvTree)
library(dplyr)

#####original extdata#####

  # Input (Example by demodata)
    file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
    Sample <- changeFormat(file = file_path)
    
  # Single call execution
    cnvTree_scDNAclustering(input = Sample,
                            pqArm_file = "hg38")
  
#####One Step#####

  setwd("/home/rstudio/yuki/")

  UCSC_cytoband_file_path <- system.file("extdata", "hg19_cytoBand.txt.gz", package = "cnvTree")

  file_path <- "ovarian_for_cnvTree/DNA/OV2295_bins_copynumber.rds"
  Sample <- changeFormat(file = file_path, cores = 5)

  cnvTree_scDNAclustering(input = Sample,
                        pqArm_file = UCSC_cytoband_file_path,
                        output.dir = "/home/rstudio/yuki/")

#####Step by Step#####
  
# Inputs
  setwd("/home/rstudio/yuki/")
  
  file_path <- "ovarian_for_cnvTree/DNA/OV2295_bins_copynumber.rds"
  Sample <- changeFormat(file = file_path, cores = 5)

# Steps for scDNA cell clustering
  pqArm_result <- pqArmClustering(input = Sample, pqArm_file = "hg19")
  Consolidation_result <- clusterConsolidation(input = Sample, pqArm_output = pqArm_result, 
                                               pqArm_file = "hg19", difratio_chr = 0.5)
  #Reclustering_output <- Reclustering(input = Sample, pqArm_output = pqArm_result,
                                    #pqArm_file = "hg38", difratio_chr = 0.5)
                                                               # default: difratio_chr = 0.3
  Subclone_output <- SubClustering(input = Sample, Consolidating_output = Consolidation_result,
                                        min_cell = 5, overlap_region = 10^8, dif_ratio = 0.5)
                             # default: min_cell = 5, overlap_region = 10^7, dif_ratio = 0.2

#combine two functions for grid search
  # combinedclustering <- function(input = Sample, pqArm_output = pqArm_result, pqArm_file = "hg38", 
  #                                difratio_chr, min_cell = 5, 
  #                                overlap_region, dif_ratio) {
  #   
  #   # Step 1: Reclustering
  #   Reclustering_output <- Reclustering(input = input, 
  #                                       pqArm_output = pqArm_output, 
  #                                       pqArm_file = pqArm_file, 
  #                                       difratio_chr = difratio_chr)
  #   
  #   # Step 2: Subclone Clustering
  #   Subclone_output <- SubcloneClustering(input = input, 
  #                                         Reclustering_output = Reclustering_output, 
  #                                         min_cell = min_cell, 
  #                                         overlap_region = overlap_region, 
  #                                         dif_ratio = dif_ratio)
  #   
  #   return(Subclone_output)  # 輸出最終結果
  # }
  

#gridSearch
  # library(NMOF)
  # gridsearch_output <- gridSearch(fun = combinedclustering, 
  #                                 levels = list(difratio_chr = seq(0, 0.5, length.out = 5),#length.out, 想產出這個範圍的多少個值
  #                                                    overlap_region = seq(10^7, 10^8, length.out = 5),
  #                                                    dif_ratio = seq(0, 0.5, length.out = 5)), 
  #                                 input = Sample, pqArm_output = pqArm_result, pqArm_file = "hg38",min_cell = 5,
  #                                 printDetail = TRUE)
  
  
# Outputs the cell clustering result
  setwd("/home/rstudio/yuki/")
  scDNA_Output(input = Sample, Summary = Subclone_output, pqArm_file = "hg19", consecutive_region = 10^4, cellcutoff = 5)
                                                                    #default: consecutive_region = 10^7, cellcutoff = 5
