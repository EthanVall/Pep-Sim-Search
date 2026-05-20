#Steps to Run Function in R Studio
1) setwd() containing the "peptide_similarity_search.R" file downloaded from this repo
2) source("peptide_similarity_search.R")
3) results <-run_peptide_search()
4) Follow prompts in console... prints colnames for convenience
5) #See the alignment of your query vs any hit
print_alignment("SIINFEKL", results$top_similarity$Peptide.Sequence[1])
print_alignment("SIINFEKL", results$top_identity$Peptide.Sequence[1])

****EXAMPLE OUTPUT****
Alignment  (| identical  + conservative  . non-conservative)
Query    : L L D A D L F F L 
         : + + . . . + . | | 
Candidate: I V V D T I M F L 

  Identity   : 2/9 (22.2%)
  Similarity : 5/9 (55.6%)  [identical + conservative]
