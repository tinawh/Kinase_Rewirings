# hypermodules prep
require(data.table); require(tidyr); require(dplyr); require(reshape2)

# location
loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/"

# variables
activedriver.cutoff <- 0.05

# files 
clinical.data <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/clinical_PANCANatlas_patient_with_followup.tsv")
longest.refseqs <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/longest_refseqs.rds")
tss <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/tss.rds")
activedriver.results <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/pval_results_all.txt")
filtered.patients <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/filtered_patients.rds")

# file names
clinical.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/clinical_data/"
unsplit.clinical.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_fixed_clin.rds"
mut.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/mutation_data/"
unsplit.muts.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_mutation_data.rds"


############################################################################################################################################
GetMutationDataByCancerType <- function(clinical.data, longest.refseqs, tss, filtered.patients, activedriver.results) {
	# Create mutation dataset for hypermodules input 
	clin <- clinical.data %>% 
				mutate(vital_status = replace(vital_status, which(vital_status %in% c("[Not Available]", "[Not Applicable]", "[Discrepancy]")), NA),
			   		days_to_last_followup = as.numeric(days_to_last_followup),
			   		days_to_death = as.numeric(days_to_death)) %>%	
				select(bcr_patient_barcode, vital_status, days_to_last_followup, days_to_death) %>%
				mutate(surv_time = pmax(days_to_last_followup, days_to_death, na.rm = T)) %>% 
				select(bcr_patient_barcode, vital_status, surv_time) %>% 
				mutate(vital_status = gsub("Dead", "DECEASED", vital_status), cancer_type = substr(bcr_patient_barcode, 6, 7)) %>% #derive cancer_type from tcga barcodes
				mutate(vital_status = gsub ("Alive", "ALIVE", vital_status)) %>%
				left_join(tss, by = c("cancer_type" = "TSS_Code")) %>% # get cancer names
				select(bcr_patient_barcode, vital_status, surv_time, Study_Name) %>% 
				mutate(Study_Name = gsub(" ", "_", Study_Name)) %>%
				right_join(filtered.patients, by = c("bcr_patient_barcode" = "Patient")) %>% # filter out hypermutated patients
				filter(!is.na(surv_time) & !is.na(vital_status) & !is.na(Study_Name) & !is.na(bcr_patient_barcode)) %>%
				filter(surv_time > 0) %>% # changed from >= 0 to > 0 
				unique()

	# create mutated patients table based on patients with available clinical data
	joined <- longest.refseqs %>%
				select(patient, hugo_name) %>%
				mutate(patient = strsplit(patient, split = ";")) %>% 
				unnest(patient) %>%
				inner_join(clin, by = c("patient" = "bcr_patient_barcode")) %>%
				select(patient, hugo_name, Study_Name) %>%
				unique()

	filtered.results <- activedriver.results %>%
							filter(fdr <= activedriver.cutoff) %>%
							select(gene) %>% 
							inner_join(joined, by = c("gene" = "hugo_name")) %>%
							na.omit() %>%
							unique()
	return(filtered.results)
}

WriteMutationDataByCancerType <- function(filtered.results, mut.folder) {
	# write mutation datasets to text files for input
	dir.create(mut.folder)

	mutation.hugo <- split(filtered.results, filtered.results$Study_Name) 

	for (dframes in names(mutation.hugo)) { 
		write.table(mutation.hugo[[dframes]][1:2], file = paste0(mut.folder, "mut_", dframes, "_hugo.csv"), quote = F, sep = ",", row.names = F, col.names = F)
}

GetClinicalDataByCancerType <- function(clinical.data, filtered.results) {
		
	clin <- clinical.data %>% 
				mutate(vital_status = replace(vital_status, which(vital_status %in% c("[Not Available]", "[Not Applicable]", "[Discrepancy]")), NA),
			   		days_to_last_followup = as.numeric(days_to_last_followup),
			   		days_to_death = as.numeric(days_to_death)) %>%	
				select(bcr_patient_barcode, vital_status, days_to_last_followup, days_to_death) %>%
				mutate(surv_time = pmax(days_to_last_followup, days_to_death, na.rm = T)) %>% 
				select(bcr_patient_barcode, vital_status, surv_time) %>% 
				mutate(vital_status = gsub("Dead", "DECEASED", vital_status), cancer_type = substr(bcr_patient_barcode, 6, 7)) %>% #derive cancer_type from tcga barcodes
				mutate(vital_status = gsub ("Alive", "ALIVE", vital_status)) %>%
				left_join(tss, by = c("cancer_type" = "TSS_Code")) %>% # get cancer names
				select(bcr_patient_barcode, vital_status, surv_time, Study_Name) %>% 
				mutate(Study_Name = gsub(" ", "_", Study_Name)) %>%
				right_join(filtered.patients, by = c("bcr_patient_barcode" = "Patient")) %>% # filter out hypermutated patients
				filter(!is.na(surv_time) & !is.na(vital_status) & !is.na(Study_Name) & !is.na(bcr_patient_barcode)) %>%
				filter(surv_time > 0) %>% # changed from >= 0 to > 0 
				unique()

	filtered.clin <- filtered.results %>% 
						select(patient) %>% 
						inner_join(clin, by = c("patient" = "bcr_patient_barcode")) %>%
						unique()

	return(filtered.clin)
}

WriteClinicalDataByCancerType <- function(filtered.clin, clinical.folder) {
	dir.create(clinical.folder)
	clin.by.cancer <- split(filtered.clin, filtered.clin$Study_Name)

	for (dframe in names(clin.by.cancer)) {
		write.table(clin.by.cancer[[dframe]][1:3], file = paste0(clinical.folder, dframe, "_clin.csv"), quote = F, sep = ",", row.names = F, col.names = F)
	}
}

#######################################################################################################################
filtered.results <- GetMutationDataByCancerType(clinical.data, longest.refseqs, tss, filtered.patients, activedriver.results)
saveRDS(filtered.results, unsplit.muts.file)

WriteMutationDataByCancerType(filtered.results, mut.folder)

filtered.clin <- GetClinicalDataByCancerType(clinical.data, filtered.results)
saveRDS(filtered.clin, unsplit.clinical.file)

WriteClinicalDataByCancerType(filtered.clin, clinical.folder)







