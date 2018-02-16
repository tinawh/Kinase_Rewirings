# Map mimp results to hypermodules for later mutation matrices 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(survival); require(survminer)

options(warn = 1)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-05-18/"
dir.create(loc)

# constants
hypermodules.dir <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/hypermodules_100/"
hypermodules.prefix <- "hypermodules_100_" # number of permutations
patient.threshrold <- 5 # threshold for number of patients in a module
coxph.threshold <- 0.05

# files
nosubset.results <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/nosubset_hypermodules_results.rds")
mimp <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/mimp.results.rds")

# file names
mimp.added.modules.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-05-18/mimp_added_modules.rds"
mimp.modules.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-05-18/mimp_modules.rds"

MapMIMPToHypermodulesPatients <- function(mimp, nosubset.results) {
	# Map hypermodules patients to mimp mutations 

	# unnest patients 
	modules <-	nosubset.results %>%
						mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>%
						unnest(List_of_patients)

	# join modules with mimp 
	mimp.joined <- mimp %>% 
						select(patient, hugo_name, mut, pwm, effect) %>% 
						inner_join(modules, by = c("patient" = "List_of_patients")) %>%
						unique()

	return(mimp.joined)
}

FindMIMPGenes <- function(mimp.joined) {
	# Match mimp genes
	
	mod.lst <- split(mimp.joined, f = mimp.joined$cancer)

	full.lsted <- list()
	for (cancer in names(mod.lst)) {

		# format cancer module
		formatted <- mod.lst[[cancer]] %>%  
						mutate(mutation = paste(hugo_name, mut)) %>% 
						select(Module, patient, hugo_name, mutation, pwm) %>% 
						unique() %>%
						mutate(module_gene = strsplit(Module, split = ";")) %>% 
						unnest(module_gene)

		gene.splits <- formatted %>%
							filter(pwm == module_gene) %>% 
							mutate(Module = strsplit(Module, split = ";"), gene_pwm = paste(hugo_name, pwm, sep = ";")) %>%
							mutate(gene_pwm = strsplit(gene_pwm, split = ";"))

			if(nrow(gene.splits) == 0) {
					full.lsted[[cancer]] <- FALSE
				} else {
					one <- gene.splits$Module
					names(one) <- as.character(seq(1:length(one)))
					eight <- gene.splits$gene_pwm
					names(eight) <- as.character(seq(1:length(eight)))
					lsted <- list()
					for (item in names(one)) {
						lsted[[item]] <- all(eight[[item]] %in% one[[item]])
					}

			if (any(lsted) == FALSE) {
				full.lsted[[cancer]] <- FALSE } else {
				full.lsted[[cancer]] <- gene.splits[which(lsted == TRUE), ]
			}
		}
	}
	return(full.lsted)	
}			

CollapseMIMPModules <- function(module.hits) {
	# Get modules with mimp mutation pairs 

	# remove cancers with no hits
	for (module in names(module.hits)) {
		if (class(module.hits[[module]]) != "data.frame") {
			module.hits[[module]] <- NULL
		}
	}

	collapsed <- list()
	for (cancer in names(module.hits)) {
		mods <- module.hits[[cancer]] %>% 
					mutate(Module = unlist(lapply(Module, paste, collapse = ";"))) %>% 
					select(Module, patient, mutation, pwm) %>% 
					mutate(patient_mut_rewiring = paste(mutation, pwm, sep = "~")) %>% 
					unique() %>% 
					group_by(Module) %>% 
					summarize(patient_muts = paste(patient_mut_rewiring, collapse = ";"), patient = paste(patient, collapse = ","), number_of_mimp_patients = n())
		collapsed[[cancer]] <- mods
 	}
 	binded <- bind_rows(collapsed, .id = "cancer")
 	return(binded)
}

AddMIMPAnnotations <- function(mimp.modules, nosubset.results) {
	# Add mimp annotations to the full set of hypermodules modules 

	mimp.selected <- mimp.modules %>% 
							select(cancer, Module, number_of_mimp_patients) %>% 
							mutate(mimp = 1)

	mod.joined <- left_join(nosubset.results, mimp.selected, by = c("cancer", "Module"))
	mod.joined[is.na(mod.joined)] <- 0

	return(mod.joined)
}
#########################################################################################################
mimp.joined <- MapMIMPToHypermodulesPatients(mimp, nosubset.results)
module.hits <- FindMIMPGenes(mimp.joined)
mimp.modules <- CollapseMIMPModules(module.hits) # 34 mimp modules out of 101
saveRDS(mimp.modules, mimp.modules.file)

mimp.added.modules <- AddMIMPAnnotations(mimp.modules, nosubset.results)
saveRDS(mimp.added.modules, mimp.added.modules.file) 

