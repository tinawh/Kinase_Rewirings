# Use biogrid network to run hypermodules with full mutations 

require(data.table); require(tidyr); require(dplyr); require(reshape2)
require(survival); require(survminer); require(stringr); require(GSA)

# constants

# location 
loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/"
dir.create(loc)

# files 
gmt <- GSA.read.gmt("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/hsapiens.pathways.NAME.gmt") 
biogrid <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/PPI_network_BioGrid_HC.txt")
clinical.data <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/clinical_PANCANatlas_patient_with_followup.tsv")
longest.refseqs <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/longest_refseqs.rds")
tss <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/tss.rds")
filtered.patients <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/filtered_patients.rds")

# file names
biogrid.network.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv"
clinical.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/clinical_data/"
unsplit.clinical.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_fixed_clin.rds"
mut.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/mutation_data/"
unsplit.muts.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_mutation_data.rds"

WriteBiogridNetwork <- function(gmt, biogrid, fil = biogrid.network.file) {
	# Combine gmt with biogrid to only include genes from gmt 

	gmt.genes <- unique(unlist(gmt$genesets))
	multinames <- biogrid %>% 
						filter(grepl("//", protein1))

	filtered <- biogrid %>% 
					mutate(protein1 = gsub (" //.*", "", protein1)) %>%
					mutate(protein2 = gsub (" //.*", "", protein2)) %>%
					filter(protein1 %in% gmt.genes & protein2 %in% gmt.genes) %>% 
					select(protein1, protein2) %>%
					unique()

	write.table(filtered, fil, quote = F, row.names = F, col.names = F, sep = "\t")
}

GetMutationDataByCancerType <- function(clinical.data, longest.refseqs, tss, filtered.patients) {
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
				na.omit() %>%
				unique()
				
	joined <- joined[, c("hugo_name", "patient", "Study_Name")]

	return(joined)
}

WriteMutationDataByCancerType <- function(filtered.results, mut.folder) {
	# write mutation datasets to text files for input
	dir.create(mut.folder)

	mutation.hugo <- split(filtered.results, filtered.results$Study_Name) 

	for (dframes in names(mutation.hugo)) { 
		write.table(mutation.hugo[[dframes]][1:2], file = paste0(mut.folder, "mut_", dframes, "_hugo.csv"), quote = F, sep = ",", row.names = F, col.names = F)
	}
}


GetClinicalDataByCancerType <- function(clinical.data, filtered.results) {
	# Create clinical dataset for hypermodules input

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
	# Write clinical datasets to text files for input

	dir.create(clinical.folder)
	clin.by.cancer <- split(filtered.clin, filtered.clin$Study_Name)

	for (dframe in names(clin.by.cancer)) {
		write.table(clin.by.cancer[[dframe]][1:3], file = paste0(clinical.folder, dframe, "_clin.csv"), quote = F, sep = ",", row.names = F, col.names = F)
	}
}

##################################################################################################
WriteBiogridNetwork(gmt, biogrid)

filtered.results <- GetMutationDataByCancerType(clinical.data, longest.refseqs, tss, filtered.patients)
saveRDS(filtered.results, unsplit.muts.file)
WriteMutationDataByCancerType(filtered.results, mut.folder)

filtered.clin <- GetClinicalDataByCancerType(clinical.data, filtered.results) 
WriteClinicalDataByCancerType(filtered.clin, clinical.folder)

