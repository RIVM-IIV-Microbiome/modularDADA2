#' Infer Sequence Variants (Big data version)
#'
#' @name inferSequenceVariantsBD
#'
#' @details A wrapper around the \code{dada2} big data tutorial for infering
#'          amplicon sequence variants. Before running this function the 
#'          user has to run \emp{filterTrimReads} from the file
#'          \emp{01_preprocessing.R}. Specify the location of filterd forward and
#'          reverse reads.  
#'
#' @param path_to_output_name A path to the folder where outputs are to be stored.
#'                            if multiple libraries are being run, better to store in 
#'                            a common head location. At this momement important to have 
#'                            backslash at the end (see example). 
#'
#' @param plot_errors Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'
#' @param path_fwd_files_filterd Path to filtered forward files to search and use
#'                         as input.  
#'
#' @param path_rev_files_filterd Path to filtered reverse files to search and use
#'                         as input.  
#'                         
#' @param nbases_for_errors Change nbases to include at least 20% of your samples.
#'                          See \code{\link[dada2]{learnErrors}}
#' 
#' @param randomize_for_errors Logical. Default is TRUE.
#'                             See \code{\link[dada2]{learnErrors}}
#'
#' @param seed_number A random seed number passed to set.seed.  
#' 
#' @param bimera_method See \code{\link[dada2]{removeBimeraDenovo}} 
#'                         
#' @param prefix Will this to output figures, table, seqtabs.  
#' 
#' @param multithread Multiple threads to use. Default is TRUE.  
#'
#' @return A list of summary table, raw seqtab and chimera removed seqtab
#' @examples 
#  path_to_output_name = "raw_data/isv_out/"
# path_fwd_files_filterd <- "raw_data/FWD/filtered"
# path_rev_files_filterd <- "raw_data/REV/filtered"
inferSequenceVariantsBD <- function(path_to_output_name=NULL,
                                    path_fwd_files_filterd = NULL,
                                    path_rev_files_filterd = NULL,
                                    seed_number = 100,
                                    plot_errors = TRUE,
                                    multithread=TRUE,
                                    nbases_for_errors = 1e9,
                                    randomize_for_errors = TRUE,
                                    bimera_method = "consensus",
                                    prefix = "Library_1") {
  require(dada2)
  require(ggplot2)
  require(readr)
  
  if(is.null(path_fwd_files_filterd) | is.null(path_rev_files_filterd) | is.null(path_to_output_name)){
    stop("Please provide file paths in path* option(s)")
  }
  
  dir.create(path_to_output_name)
  
  filtred.fwd.files <- list.files(path_fwd_files_filterd,
                                  full.names = TRUE)
  filtred.rev.files <- list.files(path_rev_files_filterd, 
                                  full.names = TRUE)
  
  message("Number of files\n", "Forward: ", length(filtred.fwd.files), "\nReverse: ", length(filtred.rev.files))
  if(length(filtred.fwd.files)!=length(filtred.rev.files)){
    stop("No. of Forward and Reverse files do not match!!!")
  }
  
  
  #################################################################################
  ## Most likely need to change if standard RIVM out is not expected             ##
  
  run <-sapply(strsplit(basename(filtred.fwd.files), "_"), `[`, 1)
  sample <- sapply(strsplit(basename(filtred.fwd.files), "_"), `[`, 2)
  sample.names <- mapply(paste0, run, sep = "_", sample)
  sample.namesR <- mapply(paste0, run, sep = "_", sample)
  
  #################################################################################
  
  if(!identical(sample.names, sample.namesR)) {
    stop("Forward and reverse files names do not match.")
  }
  names(filtred.fwd.files) <- sample.names
  names(filtred.rev.files) <- sample.namesR
  
  
  set.seed(seed_number)
  
  # Learn forward error rates
  errF <- learnErrors(filtred.fwd.files, 
                      nbases=nbases_for_errors, 
                      multithread=multithread, 
                      randomize=randomize_for_errors) #change nbases to include at least 20% of your samples
  # Learn reverse error rates
  errR <- learnErrors(filtred.rev.files, 
                      nbases=nbases_for_errors,
                      multithread=multithread, 
                      randomize=randomize_for_errors)

  
  if(plot_errors){
    # rev folder
    # save 
#    plot_folder <- paste0(path_to_output_name,"plotErrors")
 #   dir.create(plot_folder)
    # It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
    pF <- plotErrors(errF, nominalQ=TRUE) 
    file.name.plotF <- paste0(path_to_output_name, prefix,"_", "FWD_Errors_plot.pdf")
    ggsave(file.name.plotF, h=6,w=8)
    
    pR <- plotErrors(errR, nominalQ=TRUE)
    file.name.plotR <- paste0(path_to_output_name, prefix,"_",  "REV_Errors_plot.pdf")
    ggsave(file.name.plotR, h=6,w=8)
    
    message("Estimated error rates plots stored here: \n", file.name.plotF,"\n", file.name.plotR)
  }
  
  # Check if the plotted error model looks like a good fit
  fwd.convergence <- dada2:::checkConvergence(errF) #should reach 0
  rev.convergence <- dada2:::checkConvergence(errR) #should reach 0
  
  if(fwd.convergence[length(fwd.convergence)] != 0) {
    stop("checkConvergence not zero for Learn forward error rates !!!")
  }
  
  if(rev.convergence[length(rev.convergence)] != 0) {
    stop("checkConvergence not zero for Learn reverse error rates!!!")
  }
  
  
  
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  ddFs <- vector("list", length(sample.names))
  names(ddFs) <- sample.names
  ddRs <- vector("list", length(sample.names))
  names(ddRs) <- sample.names
  
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtred.fwd.files[[sam]])
    ddF <- dada(derepF, err=errF, multithread=multithread)
    ddFs[[sam]] <- ddF ###NEW: Added this for all ddFs
    derepR <- derepFastq(filtred.rev.files[[sam]])
    ddR <- dada(derepR, err=errR, multithread=multithread)
    ddRs[[sam]] <- ddR ###NEW: Added this for all ddRs
    merger <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate = FALSE)
    mergers[[sam]] <- merger
  }
  
  seqtab <- makeSequenceTable(mergers)
  
  seqtab_name <- paste0(path_to_output_name, prefix,"_", "seqtab.rds")
  saveRDS(seqtab, seqtab_name) # CHANGE ME to where you want sequence table saved
  message("Dada2 seqtab stored here: ", seqtab_name)
  
  
  # Remove Bimera/Chimera
  seqtab_nochimera <- removeBimeraDenovo(seqtab, 
                                         method=bimera_method, 
                                         multithread=multithread)
  
  seqtab_nochimera_name <- paste0(path_to_output_name, prefix,"_", "seqtab_nochimera.rds")
  saveRDS(seqtab_nochimera, seqtab_nochimera_name) # CHANGE ME to where you want sequence table saved
  message("Bimera removed seqtab stored here: ", seqtab_nochimera_name)
  
  
  
  
  GetN <- function(x) sum(getUniques(x))
  
  infer.SV.summary <- tibble::tibble(row.names = sample.names, 
                                     dadaF = sapply(ddFs, GetN),  # Now this works
                                     dadaR = sapply(ddRs, GetN),  # Now this works
                                     merged = sapply(mergers, GetN),
                                     no.chimera = rowSums(seqtab_nochimera),
                                     run.number = prefix)
  
  
  
#  SV.summary_folder <- paste0(path_to_output_name,"inferSVSummary")
 # dir.create(SV.summary_folder)
  
  file.sv.sum.table <- paste0(path_to_output_name, prefix,"_", "ASVSummary.tsv")
  
  write.table(infer.SV.summary, file.sv.sum.table, sep = "\t")
  message("Dada2 summary stored here: ", file.sv.sum.table)
  
  return(list(isvSummaryTable= infer.SV.summary,
              seqTab = seqtab,
              seqTabNonChimera = seqtab_nochimera))
  
  
}
  
