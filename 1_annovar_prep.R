# mutation data prep for MIMP analysis
require(data.table); require(tidyr); require(dplyr); require(reshape2)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/"
dir.create(loc)

# files
maf.file <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/full_maf.bed") #3600963 lines 

# file names 
annovar.input.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/annovar_input.avinput"
filtered.patients.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/filtered_patients.rds" 

FilterMafMutations <- function(maf.file) {
	# filter out everything over 9000 mutations and "PASS" mutations
	# prepped for annovar input
	filtered.patients <- maf.file %>% 
							filter(FILTER == "PASS") %>% 
							mutate(Patient = substr(Tumor_Sample_Barcode, 1, 12)) %>%
							group_by(Patient) %>% #  patient number same as tumor sample number 
							summarize(count = n()) %>%
							select(Patient, count) %>% 
							filter(count <= 9000) %>%
							select(Patient) %>% 
							unique()

	# saved for 6_hypermodules_prep							
	saveRDS(filtered.patients, filtered.patients.file)

	temp <- maf.file %>%
					mutate(Patient = substr(Tumor_Sample_Barcode, 1, 12)) %>%
					inner_join(filtered.patients, by = "Patient") %>% 
					select(Patient, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>% 
					unique() 

	# rearrange columns
	temp <- temp[, c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Patient")]

	#aggregate Patient
	filtered <- temp %>% 
					group_by(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>% 
					summarize(Patient = paste(Patient, collapse = ";")) 

	filtered$Chromosome <- as.integer(filtered$Chromosome) 

	return(filtered)
}

###########################################################################################################################################
annovar.input <- FilterMafMutations(maf.file)
write.table(annovar.input, file = annovar.input.file, sep = "\t", quote = F, col.names = F, row.names = F) 
