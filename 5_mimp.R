# prep mimp mutations, phosphosite table, sequence data, and run mimp

require(data.table); require(tidyr); require(dplyr); require(reshape2);
require(rmimp); require(GenomicRanges); require(data.table); require(Biostrings); require(parallel)

# files 
longest.refseqs <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/longest_refseqs.rds")
PTMs <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/PTM_site_table.rds")
load("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/all_RefGene_proteins.fa.rsav")


# file names 
mimp.mutation.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/mimp_mutations.txt"
ptm.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/PTM_site_table.txt"
sequence.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/sequence.fa"
mimp.results.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/mimp.results.rds"

WriteMimpMutationFile <- function(longest.refseqs, mimp.mutation.file) {
	muts <- longest.refseqs %>%
				select(longest_refs) %>%
				filter(longest_refs != "-Inf") %>%
				unique() %>%
				separate(longest_refs, c("gene", "mutation"), sep = " ", remove = TRUE)
	write.table(muts, mimp.mutation.file, row.names = F, col.names = F, quote = F, sep = "\t")
 }

 WriteMimpPhosphositeTableFile <- function(PTMs, ptm.file) {
 	filtered <- PTMs %>%
 					select(gene, position) %>%
 					unique() 
 	write.table(filtered, ptm.file, row.names = F, col.names = F, quote = F, sep = "\t")
 }

RidLastChar <- function(row) {
	substr(row, 1, nchar(row) - 1)
}

WriteMimpSequenceFile <- function(fa2, sequence.file) {
 	nam <- paste0(">", names(fa2))
	out <- c(t(cbind(nam, fa2)))
	writeLines(out, sequence.file) 
 }

RunMIMP(mimp.mutation.file, sequence.file, ptm.file) {
	mimp.results <- mimp(mimp.mutation.file, sequence.file, ptm.file, display.results = FALSE)
}

AddHugoAndPatientsMIMP <- function(mimp.results, longest.refseqs) {
	joined <- longest.refseqs %>%
				select(patient, longest_refs, hugo_name) %>%
				filter(longest_refs != "-Inf") %>%
				separate(longest_refs, c("gene", "mut"), sep = " ", remove = TRUE) %>%
				inner_join(mimp.results, by = c("gene", "mut")) 
	return(joined)
}

#################################################################################################################
WriteMimpMutationFile(longest.refseqs, mimp.mutation.file)
WriteMimpPhosphositeTableFile(PTMs, ptm.file)

fa2 <- unlist(sapply(fa1, RidLastChar))
WriteMimpSequenceFile(fa2, sequence.file)

mimp.results <- RunMIMP(mimp.mutation.file, sequence.file, ptm.file)
joined <- AddHugoAndPatientsMIMP(mimp.results, longest.refseqs)
saveRDS(joined, mimp.results.file)
