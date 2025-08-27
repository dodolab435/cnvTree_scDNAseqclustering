
PreprocessingFile <- function(input_file) {

  file_ext <- tools::file_ext(input_file)

  if (file_ext == "txt") {
    # 強制先讀成 character，避免 numeric/character 混亂
    df <- read.table(
      input_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
  } else if (file_ext == "rds") {
    rawInput <- readRDS(input_file)
    df <- rawInput[, c("cellID", "seqnames", "start", "end", "copy.number")]
  } else {
    stop("Unsupported file type: ", input_file)
  }


  df$start       <- as.integer(df$start)
  df$end         <- as.integer(df$end)
  df$copy.number <- as.integer(df$copy.number)


  if (startsWith(as.character(df$seqnames[1]), "chr") ||
      grepl("^chr", as.character(df$seqnames[1]), ignore.case = TRUE)) {

    vec <- unique(df$seqnames)
    chr_num <- vec[grepl("^chr\\d+$", vec, ignore.case = TRUE)]
    nums <- as.numeric(gsub("(?i)^chr", "", chr_num))
    max_num <- max(nums, na.rm = TRUE)
    main_chr <- paste0("chr", 1:max_num)
    extra_chr <- setdiff(vec, main_chr)

    df$seqnames <- factor(df$seqnames, levels = c(main_chr, sort(extra_chr)))
  } else {
    vec <- unique(df$seqnames)
    nums <- as.numeric(vec[grepl("^\\d+$", vec)])
    max_num <- max(nums, na.rm = TRUE)
    main_chr <- as.character(1:max_num)
    extra_chr <- setdiff(vec, main_chr)

    df$seqnames <- factor(df$seqnames, levels = c(main_chr, sort(extra_chr)))
  }


  df <- df[order(df$cellID, df$seqnames, df$start), ]


  output_file <- paste0(tools::file_path_sans_ext(basename(input_file)), ".sorted.rds")
  saveRDS(df, file = output_file)
}

################################################################################

setwd("/home/rstudio/yuki/cnvTree_GitHub/test_txt/")
input_file <- "pre-reori_chr_example_data.txt"
PreprocessingFile(input_file)

library(cnvTree)
file_path <- "pre-reori_chr_example_data.sorted.rds"

Sample <- changeFormat(file = file_path, cores = 5)

pqArm_result <- pqArmClustering(input = Sample, pqArm_file = "hg19")
Consolidation_result <- clusterConsolidation(input = Sample, pqArm_output = pqArm_result,
                                             pqArm_file = "hg19") #, difratio_chr = 0.5

Subclone_output <- SubClustering(input = Sample, Consolidating_output = Consolidation_result,
                                 min_cell = 5) #, dif_ratio = 0.5, overlap_region = 10^6

scDNA_Output(input = Sample, Summary = Subclone_output, pqArm_file = "hg19", cellcutoff = 5,
             sexchromosome_plot=TRUE) #, consecutive_region = 10^4



