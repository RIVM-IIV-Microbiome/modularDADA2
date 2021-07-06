# Biomeutilities

summarizeTaxonomicAssignments <- function(x) {
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  taxadf <- tax_table(x) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  
  whichNA <- apply(taxadf, 2, function(x) {
    sum(is.na(x))
  })
  
  
  whichNonNA <- apply(taxadf, 2, function(x) {
    sum(!is.na(x))
  })
  
  ## dada2 especially species assignments can be
  if (any(colnames(taxadf) == "Species")) {
    whichAmbiguous <- sum(grepl("/", taxadf$Species))
    message(paste0("In Species column- ", whichAmbiguous, " are ambigous assignments"))
  }
  
  
  total <- whichNA + whichNonNA
  
  TaxonomicAssignments <- tibble::tibble(
    Ranks = rank_names(x),
    TaxonomyAssigned = whichNonNA,
    Total = total,
    PrecentAssigned = round(100 * whichNonNA / total, 1)
  )
  return(TaxonomicAssignments)
}



##############################################################################################################
#' Build a phylogenetic tree
#'
#' @name buildSeqTree
#'
#' @details A wrapper to build phylogenetic tree for ASV sequences
#'          obtained after \code{dada2}.
#'
#' @param x Sequences DNAStringSet
#'
#' @param pml_model See optim.pml \code{\link[phangorn]{pml.control}}
#'
#' @param pml_optInv See optim.pml \code{\link[phangorn]{pml.control}}
#'
#' @param pml_optGamma See optim.pml \code{\link[phangorn]{pml.control}}
#'
#' @param pml_rearrangement See optim.pml \code{\link[phangorn]{pml.control}}
#' #'
#' @param pml_control See optim.pml \code{\link[phangorn]{pml.control}}
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#' # select few for example
#' seqs <- sample(SprockettTHData@refseq, 5)
#' sp.tree <- buildSeqTree(seqs)
#' sp.tree$tree
#'
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings DNAStringSet
#' @importFrom phangorn phyDat dist.ml NJ pml optim.pml pml.control
#' @importFrom stats update
#'
#' @export

buildSeqTree <- function(x,
                         pml_model="GTR",
                         pml_optInv=TRUE,
                         pml_optGamma=TRUE,
                         pml_rearrangement = "stochastic",
                         pml_control = pml.control(trace = 0)) {
  require(phangorn)
  # global vars
  alignment <- phangAlign <- dm <- treeNJ <- fit <- fitGTR <- NULL
  
  if(class(x)=="DNAStringSet"){
    x <- x
  } else {
    x <- DNAStringSet(x)
  }
  
  alignment <- DECIPHER::AlignSeqs(x, anchor=NA, verbose=T)
  
  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
  
  dm <- phangorn::dist.ml(phangAlign)
  
  # Note, tip order != sequence order
  treeNJ <- phangorn::NJ(dm)
  
  fit <- phangorn::pml(treeNJ, data=phangAlign)
  
  fitGTR <- stats::update(fit, k=4, inv=0.2)
  
  fitGTR <- phangorn::optim.pml(fitGTR,
                                model = pml_model,
                                optInv = pml_optInv,
                                optGamma = pml_optGamma,
                                rearrangement = pml_rearrangement,
                                control = pml_control)
  
  return(fitGTR)
  
}


##############################################################################################################

#' Calculate QC Metrics for Taxa and Samples
#'
#' @name calculateQC
#'
#' @details Calculate QC metrics for taxa and samples.
#'
#' @param x A phyloseq object
#'
#' @return A list with two tibbles.
#'
#' \emph{SampleQC}
#' \itemize{
#' \item{sample_total_taxa:}{Total taxa in a sample}
#' \item{sample_total_reads:}{Total reads in a sample}
#' \item{number_taxa_account_fifty_percent:}{Number of taxa account for 50 percent}
#' }
#' \emph{TaxaQC}
#' \itemize{
#' \item{taxa_counts:}{Total counts of taxa in all samples}
#' \item{taxa_mean_counts:}{Mean counts of taxa in all samples}
#' \item{taxa_dectected_samples:}{Number of samples in which taxa detected}
#' \item{taxa_prevalence:}{Percent of samples in which taxa detected}
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' qc_Data <- calculateQC(FuentesIliGutData)
#' qc_Data
#'
#' @importFrom microbiome prevalence abundances
#' @importFrom phyloseq nsamples ntaxa sample_sums taxa_sums taxa_sums
#' @importFrom phyloseq merge_phyloseq phy_tree taxa_names<-
#' @importFrom tibble tibble
#' @importFrom stats sd
#'
#' @export

calculateQC <-  function(x) {
  
  # gloabl vars
  group <- merge_phyloseq <- taxa <-taxa_counts <- taxa_mean_counts <- NULL
  taxa_stdev_counts <- NULL
  
  ## Must be a phyloseq object
  if ( !is(x, "phyloseq") ){
    stop("input must be an phyloseq object.")
  }
  
  ## the input must have few samples
  if ( nsamples(x) < 1 ){
    stop("input must have at least one sample")
  }
  
  if ( ntaxa(x) < 1 ) {
    stop("input must have at least one taxa")
  }
  
  ## See what versions of the expression data are available in the object
  if(any(sample_sums(x)==1) | any(.check_decimal(sample_sums(x)))){
    stop("Data must be counts")
  }
  
  # Create a tibbles with QC information
  qcSample <- tibble::tibble(Samples=sample_names(x),
                             # number of taxa per sample
                             sample_total_taxa = colSums(microbiome::abundances(x) != 0),
                             # Total reads/sample
                             sample_total_reads = sample_sums(x),
                             number_taxa_account_fifty_percent = microbiome::coverage(x))
  
  qcTaxa <- tibble::tibble(taxa = taxa_names(x),
                           # number of taxa per sample
                           taxa_counts = taxa_sums(x),
                           taxa_mean_counts = rowMeans(microbiome::abundances(x)),
                           taxa_stdev_counts = apply(microbiome::abundances(x), 1, sd),
                           taxa_cv = taxa_stdev_counts/taxa_mean_counts,
                           taxa_rank = rank(rowMeans(microbiome::abundances(x))),
                           taxa_dectected_samples = rowSums(microbiome::abundances(x) != 0),
                           taxa_prevalence = microbiome::prevalence(x) *100,
                           percent_of_total = taxa_counts/sum(sample_sums(x)) *100)
  
  return(list(SampleQC=qcSample, TaxaQC=qcTaxa))
  
  
}


#' @keywords internal
#' check decimal
#' \url{https://www.rdocumentation.org/packages/schoolmath/versions/0.4/topics/is.decimal}
.check_decimal <- function(x){
  
  start <- 1
  end <- length(x)+1
  while(start<end){
    y <- x[start]
    
    test <- floor(y)
    if(y==test){
      if(start==1){
        result=FALSE
      }else{
        result<- c(result,FALSE)
      }
      
    }else{
      if(start==1){
        result=TRUE
      }else{
        result <- c(result,TRUE)
      }
    }
    start <- start+1
  }
  
  return(result)
}


######################################################################################################
#' Check Trend in Abundance and Taxonomic Assignment
#'
#' @name trendAbundanceAssignment
#'
#' @details Check if the more abundance-prevalent taxa
#'          have better taxonomic assignments. This is
#'          a pre-check, not corrected for differences
#'          in sequencing depth.
#'
#' @param x A phyloseq object
#'
#' @param quantiles Abundances values to sort. Can be changed to specify
#'                  how many values in a distribution are above or below
#'                  a certain limit.
#'
#' @param plot Logical. Default is TRUE.
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#'
#' p1 <- trendAbundanceAssignment(SprockettTHData,
#'                               quantiles = seq(0, 95, by = 5),
#'                               plot=TRUE)
#' p1 + ggplot2::scale_colour_brewer("", palette = "Spectral") +
#'     ggplot2::theme_minimal() +
#'     ggplot2::theme(axis.text.x = element_text(angle=90, vjust=0.5))
#'
#' @return Either a list with data.frame and plot or just
#'         ggplot object.
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Thorsten Brach, 2018. MicrobiomeX2 Pipeline.
#'        \emph{GitHub}.
#'        \url{https://github.com/TBrach/MicrobiomeX2}
#'        \url{https://github.com/TBrach/Dada_Pipel/blob/master/Generalized_Phyloseq_Analysis_New.Rmd}
#' }
#' }
#' @importFrom stats quantile
#' @importFrom phyloseq taxa_sums 'sample_data<-'
#' @importFrom data.table as.data.table melt
#' @importFrom tidyr separate
#' @importFrom ggplot2 ggplot geom_point scale_x_discrete ylab xlab theme element_text
#' @export

trendAbundanceAssignment <- function(x, quantiles = seq(0, 95, by = 10),
                                     plot=TRUE){
  
  #globalvars
  Quant <- Var1 <- Var2 <- variable <- aes <- value <- NULL
  taxadf <- tax_table(x) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  taxa.counts <- phyloseq::taxa_sums(x)
  
  abQuantiles <- stats::quantile(taxa.counts, probs = quantiles/100)
  
  which_taxa <- lapply(abQuantiles,
                       function(quant) {taxa.counts >= quant})
  
  numb_taxa <- sapply(which_taxa, sum)
  
  filtered.taxas <- lapply(which_taxa, function(indexes){
    taxadf[indexes,]
  })
  
  assignment_distributions <- lapply(filtered.taxas, .summarizeTaxonomicAssignmentsDF)
  
  pc_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
  
  rownames(pc_assigned) <- colnames(taxadf)
  
  colnames(pc_assigned) <- paste("Ab_", quantiles, "_", round(abQuantiles), "_", numb_taxa, sep = "")
  pc_assigned <- as.data.frame(pc_assigned)
  pc_assigned$Var1 <- rownames(pc_assigned)
  pc_assigned <- data.table::as.data.table(pc_assigned)
  
  pc_assigned_ldf <- data.table::melt(pc_assigned,
                                      id.vars = c("Var1"),
                                      value.name = "value") %>%
    tidyr::separate(variable, c("Type", "Quant", "abQuant", "numb_taxa"), sep="_")
  
  pc_assigned_ldf$Var1 <- factor(pc_assigned_ldf$Var1, levels = colnames(taxadf), ordered = TRUE)
  
  pc_assigned_ldf$Quant <- factor(pc_assigned_ldf$Quant, levels = as.character(quantiles), ordered = TRUE)
  
  trendplot <- ggplot2::ggplot(pc_assigned_ldf, aes(x = Quant, y = value))
  trendplot <- trendplot +
    ggplot2::geom_point(size = 2, aes(fill = Var1), 
                        color="grey50", shape=21, 
                        alpha=0.5) +
    ggplot2::scale_x_discrete(breaks = quantiles, labels = paste(quantiles, " (", numb_taxa, ")", sep = "")) +
    ggplot2::ylab("Percentage of taxa assigned") +
    ggplot2::xlab("Total counts quantile (No. of remaining taxa)") +
    ggplot2::theme(axis.text.x = element_text(angle=90, vjust=0.5))
  
  if(plot){
    return(trendplot)
  }
  return(list(pc_assigned_ldf, trendplot))
  
}


## Helper

.summarizeTaxonomicAssignmentsDF <- function(taxa){
  
  # taxa <- as.data.frame(unclass(tax_table(physeq)))
  
  countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
  countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
  # ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
  total <- countNA + countNonNA
  assignment_distribution <- data.frame(assigned = countNonNA,
                                        total = total,
                                        PC_assigned = round(100*countNonNA/total, 1))
  
}






