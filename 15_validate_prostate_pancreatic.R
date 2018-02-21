# External validation of pancreatic and prostate cancer patients with collaborators 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(stringr)

# constants

# files 
tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/modules_added_nonlinkers.rds")
tcga.clins <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_fixed_clin.rds")
tcga.muts <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_mutation_data.rds")
mimp.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-05-18/mimp_modules.rds")

# file names
pancreatic.modules.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/pancreatic_module_genes.txt"
prostate.modules.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/prostate_module_genes.txt"

WriteNonLinkerModuleGenes <- function(tcga.modules, cancer.type, file.name) {
	# Write nonlinker mutated genes in modules based on cancer type
	nonlinker.genes <- tcga.modules %>% 
							filter(cancer == cancer.type) %>% 
							ungroup() %>%
							select(intersect) %>% 
							mutate(intersect = strsplit(intersect, split = ";")) %>% 
							unnest() %>%
							unlist() %>% 
							unique()

	writeLines(nonlinker.genes, con = file.name, sep = "\t")
}

######################################################################################################
WriteNonLinkerModuleGenes(tcga.modules, cancer.type = "Pancreatic_adenocarcinoma", file.name = pancreatic.modules.file)
WriteNonLinkerModuleGenes(tcga.modules, cancer.type = "Prostate_adenocarcinoma", file.name = prostate.modules.file)
