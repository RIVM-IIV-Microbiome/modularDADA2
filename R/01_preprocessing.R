#' Filter and Trim Reads  
#'
#' @name filterTrimReads
#'
#' @details A wrapper around the \code{dada2} big data tutorial for preprocessing
#'          fastq files. This function will read files from directory specified 
#'          by the user, created new sub directory to store raw forward and reverse
#'          reads and filtered reads output of \code{\link[dada2]{filterAndTrim}}.
#'          Optionally, an aggregated quality profile will be generated using
#'          \code{\link[dada2]{plotQualityProfile}}. All the processing values for 
#'          \code{\link[dada2]{filterAndTrim}} can be specified by the user.
#'
#' @param path_to_all_files A path to the folder where raw fastq files are stored.
#'                          At this momement important to have backslash at the end
#'                          (see example). It is assumed that the primer sequneces 
#'                          are removed before running this step.  
#'
#' @param plot_qc_profile Wether to plot quality profiles. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotQualityProfile}} function.
#'
#' @param fwd_file_pattern A pattern that specifies fowward files to search and use
#'                         as input.  
#'
#' @param rev_file_pattern A pattern that specifies reverse files to search and use
#'                         as input.  
#'
#' @param flag.library.size The minimal library size that is required for not being 
#'                          highlighted in the trim summary output table.
#'                         
#' @param prefix Will this to output figures and table.  
#'
#'
#' @examples 
#' 
#' filterTrimReads(path_to_all_files = "raw_data/",
#'                 plot_qc_profile = TRUE,
#'                 fwd_file_pattern = "R1_001.fastq",
#'                 rev_file_pattern= "R2_001.fastq",
#'                 flag.library.size = 1000,
#'                 prefix = "Library_1",
#'                 truncLen=c(150,100), 
#'                 maxEE=c(2,2), 
#'                 truncQ=2, 
#'                 maxN=0, 
#'                 rm.phix=TRUE,
#'                 compress=TRUE, 
#'                 verbose=TRUE, 
#'                 multithread=TRUE)
#' 
#' 
#' importFrom dada2 plotQualityProfile filterAndTrim 
#' importFrom tibble rownames_to_column 
#' importFrom dplyr mutate
#' 
filterTrimReads <- function(path_to_all_files = NULL,
                            plot_qc_profile = TRUE,
                            fwd_file_pattern = "R1_001.fastq",
                            rev_file_pattern= "R2_001.fastq",
                            flag.library.size = 1000,
                            prefix = NULL,
                            ...) {
  
  require(dada2)
  require(ggplot2)
  require(readr)
  library(dplyr)
  
  main.folder <- path_to_all_files
  
  # fwd folder
  dir.create(paste0(main.folder, "/", "FWD"))
  # rev folder
  dir.create(paste0(main.folder, "/", "REV"))
  message("Subdirectories created here: \n", paste0(main.folder, "/", "REV") ,"\n", paste0(main.folder, "/", "FWD"))
  

  # copy fwd files to fwd folder
  fwd.folder <- paste0(main.folder, "/", "FWD")
  list.files.fwd <- list.files(main.folder,fwd_file_pattern)
  list.files.fwd <- paste0(main.folder, list.files.fwd)
  file.copy(list.files.fwd, fwd.folder)
  
  # copy REV files to REV folder
  rev.folder <- paste0(main.folder, "/", "REV")
  list.files.rev <- list.files(main.folder, rev_file_pattern)
  list.files.rev <- paste0(main.folder, list.files.rev)
  file.copy(list.files.rev, rev.folder)
  
  # get a list of fdw and rev files

  fwd.files <- sort(list.files(fwd.folder, pattern="fastq.gz"))
  rev.files <- sort(list.files(rev.folder, pattern="fastq.gz"))
  
  ## Check if pairs exists
  if(length(fwd.files) != length(rev.files)) {
    stop("Forward and reverse files do not match.")
  }
  
  if(plot_qc_profile){
    # rev folder
    # save 
    plot_folder <- paste0(main.folder,"plotQualityProfile")
    dir.create(plot_folder)
    
    # check quality profiles
    plotQualityProfile(paste0(fwd.folder,"/",fwd.files), 
                       n = 5e+05, 
                       aggregate = TRUE)
    file.name.plotF <- paste0(plot_folder, "/", prefix,"_", "FWD_QC_aggregated_plot.pdf")
    ggsave(file.name.plotF, h=6,w=8)
    
    plotQualityProfile(paste0(rev.folder,"/",rev.files), 
                       n = 5e+05, 
                       aggregate = TRUE)
    file.name.plotR <- paste0(plot_folder, "/", prefix,"_", "REV_QC_aggregated_plot.pdf")
    ggsave(file.name.plotR, h=6,w=8)
    message("QualityProfile stored here: \n", file.name.plotF,"\n", file.name.plotR)
  }
  
  # location for filtered files
  
  # create path for filtered reads FWD and REV
  fwd.files.filterd <- file.path(fwd.folder, "filtered") 
  rev.files.filterd <- file.path(rev.folder, "filtered") 
  
  # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
  trim.info <-filterAndTrim(fwd=file.path(fwd.folder, fwd.files), 
                            filt=file.path(fwd.files.filterd, fwd.files),
                            rev=file.path(rev.folder, rev.files), 
                            filt.rev=file.path(rev.files.filterd, rev.files),
                            ...)
  
  message("Filtered fastq (fastq.gz) files stored here: \n", fwd.files.filterd,"\n", rev.files.filterd)
  
  # filterAndTrim QC cals
  trim.summary <- trim.info %>% 
    as.data.frame() %>%  
    tibble::rownames_to_column("input_file") %>% 
    dplyr::mutate(reads.lost.number= reads.in - reads.out,
                  reads.lost.percent= round(((reads.in - reads.out)/reads.in *100), 2),
                  is.low.reads.in = ifelse(reads.in < flag.library.size, "YES", "NO"),
                  is.low.reads.out = ifelse(reads.in < flag.library.size, "YES", "NO"),
                  run.number = prefix)
  
  trim_folder <- paste0(main.folder,"trimSummary")
  dir.create(trim_folder)
  
  file.name.trim.table <- paste0(trim_folder, "/", prefix,"_", "trimSummary.tsv")
  
  write.table(trim.summary, file.name.trim.table, sep = "\t")
  message("Filter trim summary stored here: ", file.name.trim.table)
  
  return(trim.summary)
  
}



