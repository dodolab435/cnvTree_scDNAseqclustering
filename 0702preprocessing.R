
PreprocessingFile <- function(input_file) {

  ################################################################################

  file_ext <- tools::file_ext(input_file)

  ################################################################################

  if (file_ext == "txt"){
    noFactorInput <- read.delim2(input_file, sep = "\t", header = TRUE)
  } else if (file_ext == "rds"){
    rawInput <- readRDS(input_file)
    noFactorInput <- rawInput[, c("cellID", "seqnames", "start", "end", "copy.number")]
  } else {
    stop("Unsupported file type: ", input_file)
  }

  ################################################################################

  noFactorInput$start       <- as.integer(noFactorInput$start)
  noFactorInput$end         <- as.integer(noFactorInput$end)
  noFactorInput$copy.number <- as.integer(noFactorInput$copy.number)

  ################################################################################

  if (startsWith(as.character(noFactorInput$seqnames[1]), "chr")||
    grepl("^chr", as.character(noFactorInput$seqnames[1]), ignore.case = TRUE)){ #忽略chr大小寫

    vec <- unique(noFactorInput$seqnames)
    chr_num <- vec[grepl("^chr\\d+$", vec, ignore.case = TRUE)] #chr+數字
    nums <- as.numeric(gsub("(?i)^chr", "", chr_num)) #移除chr只留數字

    max_num <- max(nums, na.rm = TRUE)
    main_chr <- paste0("chr", 1:max_num)
    extra_chr <- setdiff(vec, main_chr) #抓非數字chr

    noFactorInput$seqnames <- factor(
      noFactorInput$seqnames,
      levels = c(main_chr, sort(extra_chr))
    )
  } else {

    vec <- unique(noFactorInput$seqnames)

    nums <- as.numeric(vec[grepl("^\\d+$", vec)])  # 純數字
    max_num <- max(nums, na.rm = TRUE)
    main_chr <- as.character(1:max_num)
    extra_chr <- setdiff(vec, main_chr)

    noFactorInput$seqnames <- factor(
      noFactorInput$seqnames,
      levels = c(main_chr, sort(extra_chr))
    )
  }

  ################################################################################

  noFactorInput <- noFactorInput[order(noFactorInput$cellID,
                                       noFactorInput$seqnames,
                                       noFactorInput$start), ]

  ################################################################################

  output_file <- paste0(tools::file_path_sans_ext(basename(input_file)), ".sorted.rds")
  saveRDS(noFactorInput, file = output_file)

}

  ################################################################################

  setwd("/home/rstudio/yuki/cnvTree_GitHub/test_txt/")

  df <- read.table(
    "reori_chr_example_data.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character"   # 每一欄都強制讀成字串
  )
  df$start <- as.numeric(df$start)
  df$end   <- as.numeric(df$end)

  write.table(df, "pre-reori_chr_example_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)




  input_file <- "pre-reori_chr_example_data.txt"

  input_file <- "txt_files/fixedOvarian.txt" # txt
  input_file <- "rds_files/fixedOvarian.rds" # rds
  PreprocessingFile(input_file)

###cnvTree original codes###

library(cnvTree)
file_path <- "pre-reori_chr_example_data.sorted.rds"
file_path <- "fixedOvarian.sorted.rds"

Sample <- changeFormat(file = file_path, cores = 5)

pqArm_result <- pqArmClustering(input = Sample, pqArm_file = "hg19")
Consolidation_result <- clusterConsolidation(input = Sample, pqArm_output = pqArm_result,
                                             pqArm_file = "hg19") #, difratio_chr = 0.5

Subclone_output <- SubClustering(input = Sample, Consolidating_output = Consolidation_result,
                                 min_cell = 5) #, dif_ratio = 0.5, overlap_region = 10^6

scDNA_Output(input = Sample, Summary = Subclone_output, pqArm_file = "hg19", cellcutoff = 5,
             sexchromosome_plot=TRUE) #, consecutive_region = 10^4



################################################################################

  # Input (Example using demo data)
  file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
  Example_data <- changeFormat(file = file_path)

  # Single-call execution
  cnvTree_scDNAclustering(input = Example_data,
                        pqArm_file = "hg38")








