
# ModularDADA2  

Converting the dada2 big data workflow into a modular pipeline for analysing multiple sequencing libraries.  

Divided into 

  * Preprocessing  
  * Infer Sequence Variants and remove chimera  
  * Taxonomic assignments

## Codes  
This is accomplished by set of three functions.  

  * Preprocessing  
      * `filterTrimReads` *(R/01_preprocessing.R)*  
      
  * Infer Sequence Variants and remove chimera   
      * `inferSequenceVariantsBD` *(R/02_inferSequenceVariantsBD.R)*  
      
  * Taxonomic assignments   
      * `classifyASV` *(R/03_taxonomyAssignments.R)*   
      
