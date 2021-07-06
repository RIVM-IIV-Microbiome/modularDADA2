#' Classify ASV 
#' 
#' @name classifyASV
#'
#' @details A wrapper around the \code{dada2} big data tutorial for classifying
#'          infered amplicon sequence variants. Before running this function the 
#'          user has to run \emp{filterTrimReads} from the file
#'          \emp{01_preprocessing.R}, \emp{inferSequenceVariantsBD} from 
#'          \emp{02_inferSequenceVariantsBD.R}. This functions takes in a seqtab either per
#'          run or from multiple runs that are merged by user 
#'          using \code{\link[dada2]{mergeSequenceTables}}.
#'
#' @param path_to_output_name A path to the folder where outputs are to be stored.
#'                            if multiple libraries are being run, better to store in 
#'                            a common head location. At this momement important to have 
#'                            backslash at the end (see example). 
#'
#' @param training_set Path to the training set to be passed on to 
#'                        calls the \code{\link[dada2]{assignTaxonomy}} function.
#'                        
#' @param training_set_species Path to the species training set to be passed on to  
#'                        calls the \code{\link[dada2]{addSpecies}} function.
#'                        
#' @param minBoot Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'                                                
#' @param multithread Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'                        
#' @param tryRC Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'                        
#' @param verbose Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'   
#' @param verbose Wether to plot estimated errors. Default is TRUE which 
#'                        calls the \code{\link[dada2]{plotErrors}} function.
#'                                                                            
#'
classifyASV <- function(input_seqtab = input_seqtab, 
                        path_to_output_name = "myoutput/",
                        training_set = "silvaDBv138.1/silva_nr99_v138.1_train_set.fa.gz",
                        training_set_species = "silvaDBv138.1/silva_species_assignment_v138.1.fa.gz",
                        minBoot=80,
                        multithread = 4,
                        tryRC = TRUE,
                        verbose = TRUE,
                        set.seed=356,
                        allowMultiple=3){
  
  require(dada2)
  require(Biostrings)
  require(tidyverse)
  
  set.seed(set.seed)
  
  if(ncol(input_seqtab) < 3000){
    #message(paste0("Dataset is too big to be processed at once. Divind dataset in ", nr_chunks, " chunks to avoid memory crashing..."))
    dna <- DNAStringSet(getSequences(input_seqtab)) 
    taxa <- assignTaxonomy(dna, 
                           refFasta = training_set, 
                           minBoot=minBoot,
                           multithread = multithread,
                           tryRC = tryRC,
                           verbose = verbose)
    
    
    taxaSpecies <- addSpecies(taxa, 
                              refFasta = training_set_species,
                              allowMultiple =allowMultiple, 
                              verbose=TRUE)
    
    tx_path <- paste0(path_to_output_name, "taxa_assignTaxonomy.rds")
    saveRDS(taxa, tx_path)
    
    tx_path2 <- paste0(path_to_output_name, "taxa_addSpecies.rds")
    saveRDS(taxaSpecies, tx_path2)
    
  } else if (ncol(input_seqtab) > 3000){
    
    # section from https://gitlab.rivm.nl/hernanda/triumph_dada2_run_pipeline/blob/master/bin/3_inferSAVs.R
    
    # Number of chunks to divide seqtab in 3000 ASVs max
    chunk_size <- 3000
    nr_chunks <- ceiling(ncol(input_seqtab)/chunk_size) 
    
    message(paste0("Dataset is too big to be processed at once. Divind dataset in ", nr_chunks, " chunks to avoid memory crashing..."))
    
    seqtab_chunks <- vector(mode = "list", length = nr_chunks)
    last_number <- 0
    for(i in 1:nr_chunks){
      
      columns_to_add <- seq(1+last_number, ifelse(chunk_size*i > ncol(input_seqtab), ncol(input_seqtab), chunk_size*i))
      seqtab_chunks[[i]] <- input_seqtab[ , columns_to_add]
      last_number <- chunk_size*i
      
      
      #Assign taxonomy per chunk
      message("Assigning taxonomy to each chunk...")
      taxa_run <- seqtab_chunks %>%
        map(assignTaxonomy, refFasta = training_set, minBoot=minBoot,
            multithread = multithread,
            tryRC = tryRC,
            verbose = verbose)
      message("Saving taxonomy RDS...")
      tx_path <- paste0(path_to_output_name, "taxa_assignTaxonomy.rds")
      saveRDS(taxa_run, tx_path)
      
      #saveRDS(taxa_run, taxonomy_file)
      
      #Adding species per chunk
      message("Adding species...")
      
      species_run <- taxa_run %>%
        map(addSpecies, refFasta = training_set_species, allowMultiple = allowMultiple) %>%
        map(as.data.frame)%>%
        bind_rows()
      rownames(species_run) <- colnames(input_seqtab)
      
      message("Saving species RDS...")
      
      tx_path2 <- paste0(path_to_output_name, "taxa_addSpecies.rds")
      saveRDS(species_run, tx_path2)
      
      
      #saveRDS(species_run, species_file)
      message("Finished successfully!")
      
    }
    
  }
 
  
  
}