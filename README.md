Simple immunopeptidomics data tool intended for finding cross reactive peptides within a reference sample
Requires a Peptide Sequence column containing letter symbols of amino acids i.e. SIINFEKL
Requires a Normal or non-diseased sample column with intensities (Any peptide with intensity > 0 will be Aligned).
With the emergence of computationally designed TCR mimentics and minibinders targeting specific pMHCs, this tool was intended to identify 
potential cross reactive peptides that posess analogous biochemical properties to a disease specific target.
Can be used without in house immunopeptidomics data with the HLA Ligand atlas. CNS specific class 1 ligand atlas data can be found in repository
Use "Peptide.Length" as input for the HLA Ligand atlas when prompted to enter a normal intensity column name. 
Query peptide does not need to be found in immunopeptidomics data set

#Steps to Run Function in R Studio
1) setwd() containing the "peptide_similarity_search.R" file downloaded from this repo
2) source("peptide_similarity_search.R")
3) results <-run_peptide_search()
4) Follow prompts in console... prints colnames for convenience
5) #See the alignment of your query vs any hit
print_alignment("SIINFEKL", results$top_similarity$Peptide.Sequence[1])
print_alignment("SIINFEKL", results$top_identity$Peptide.Sequence[1])
<img width="500" height="213" alt="image" src="https://github.com/user-attachments/assets/6a2a7698-ce7e-42b2-9892-85269bdbcd61" />


