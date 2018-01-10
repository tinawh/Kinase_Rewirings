# mutation data prep for MIMP 

require(data.table); require(tidyr); require(dplyr); require(reshape2); 
require(hash); require(magrittr); require(ggpubr); require(ggsci); 
require(gProfileR); require(gridExtra); require(cowplot);
require(rvest); require(xml2)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/"
maf.file <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/full_maf.bed") #3600963 lines 

#file names 
annovar.input.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/annovar_input.avinput")

FilterMafMutations <- function(maf.file) {
	#filter out everything over 900 mutations and "PASS" mutations
	# prepped for annovar input
	maf.file %>% 
		filter(FILTER == "PASS") %>% 
		mutate(Patient = substr(Tumor_Sample_Barcode, 1, 12)) %>%
		group_by(Patient) %>% #  patient number same as tumor sample number 
		summarize(count = n()) %>%
		select(Patient, count) %>% 
		filter(count <= 9000) %>%
		select(Patient) %>% 
		unique() -> filtered.patients #9060 patients left 

	maf.file %>%
		mutate(Patient = substr(Tumor_Sample_Barcode, 1, 12)) %>%
		inner_join(filtered.patients, by = "Patient") %>% 
		select(Patient, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>% 
		unique() -> temp 

	temp <- temp[, c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Patient")]

	temp %>% 
		group_by(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>% 
		summarize(Patient = paste(Patient, collapse = ";")) -> filtered  
	filtered$Chromosome <- as.integer(filtered$Chromosome) 
	return(filtered)
}

###########################################################################################################################################
annovar_input <- FilterMafMutations(maf.file)
write.table(annovar_input, file = annovar.input.file, sep = "\t", quote = F, col.names = F, row.names = F) 
