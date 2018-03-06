# Map mimp results to hypermodules for later mutation matrices 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(survival); require(survminer); require(stringr)

options(warn = 1)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/"
#dir.create(loc)
network.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/cytoscape_networks/" 
dir.create(network.folder)

# constants

# files
mimp.added.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/mimp_added_modules.rds")
mimp <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/mimp.results.rds")
clins <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_fixed_clin.rds")
muts <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_mutation_data.rds")
mimp.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/mimp_modules.rds")

network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)

# file names
survival.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/scurves.pdf" 
mutation.matrix.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/mutation.matrix.pdf" 
modules.non.linker.added.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/modules_added_nonlinkers.rds" 

JoinModuleMutations <- function(modules, muts) {
	# Join modules with mutation data 

	unnested.patients <- modules %>% 
								mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>%
								unnest(List_of_patients) 

	mut.joined <- inner_join(unnested.patients, muts, by = c("List_of_patients" = "patient", "cancer" = "Study_Name"))

	return(mut.joined)
}

AddNonLinker <- function(mut.joined, mimp) {
	# Add column with only non-linkers module genes 

	aggregated <- mut.joined %>% 
						group_by(Module, cancer) %>%
						summarize(List_of_patients = unique(paste(List_of_patients, collapse = ",")), gene = unique(paste(hugo_name, collapse = ";")))

	# get intersection of module and patient genes to remove linker genes 
	compared <- mapply(FUN = function(mod.genes, patient.muts) {
					intersect(mod.genes, patient.muts)
					}, mod.genes = strsplit(aggregated$Module, split = ";"), patient.muts = strsplit(aggregated$gene, split = ";"))

	# add intersect
	aggregated$intersect <- unlist(lapply(compared, function(element) paste(element, collapse = ";"))) 

	# rejoin patients and mutations
	intersected.module.genes <- select(aggregated, Module, cancer, intersect)
	
	rejoined <- inner_join(intersected.module.genes, mimp.added.modules, by = c("Module" = "Module", "cancer" = "cancer"))

	return(rejoined)
}

DrawSurvivalCurvesByCancerType <- function(nonlinkers.added, muts, clins, survival.file) {
	# Draw survival curves 

	# prep modules
	filtered <- nonlinkers.added %>% 
					select(Module, cancer, intersect, List_of_patients, coxph, logrank) 

	# split modules by cancer
	filtered.muts <- split(filtered, f = filtered$cancer)

	# split clin by cancer and extract only matching cancers
	split.clins.full <- split(clins, f = clins$Study_Name)
	split.clins <- split.clins.full[names(filtered.muts)]

	lst <- list()
	for (cancer in names(filtered.muts)) { 
		pvalue.arranged <- filtered.muts[[cancer]] %>% 
								select(Module, intersect, coxph, logrank) %>% 
								arrange(coxph) %>%
								unique() 

		filtered.muts[[cancer]]$Module <- factor(filtered.muts[[cancer]]$Module, levels = pvalue.arranged$Module)
		filtered.muts.cancer <- split(filtered.muts[[cancer]], f = filtered.muts[[cancer]]$Module) 

		for (dframe in names(filtered.muts.cancer)) {
			filtered.muts.cancer[[dframe]]$Module <- as.character(filtered.muts.cancer[[dframe]]$Module)
			Module <- unique(filtered.muts.cancer[[dframe]]$Module)
			intersect.genes <- unique(filtered.muts.cancer[[dframe]]$intersect)
			coxph.value <- formatC(unique(filtered.muts.cancer[[dframe]]$coxph), format = "e", digits = 2)
			logrank.value <- formatC(unique(filtered.muts.cancer[[dframe]]$logrank), format = "e", digits = 2)

			df.update <- filtered.muts.cancer[[dframe]] %>%
							mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
							select(Module, List_of_patients) %>%
							unnest(List_of_patients) %>%
							unique() %>%
							mutate(is_module_mutated = 1) %>%
							full_join(split.clins[[cancer]], by = c("List_of_patients" = "patient")) %>%
							mutate(vital_status = gsub("DECEASED", 1, vital_status)) %>%
							mutate(vital_status = gsub("ALIVE", 0, vital_status)) %>%
							mutate(vital_status = as.numeric(vital_status)) %>% 
							ungroup() %>%
							select(List_of_patients, is_module_mutated, vital_status, surv_time)

				df.update$is_module_mutated[is.na(df.update$is_module_mutated)] <- 0

				caption <- paste(paste("logrank p = ",logrank.value, "; coxph p = ", coxph.value, sep = ""), intersect.genes, cancer, sep = "\n")

				logrank.fit <- survfit(Surv(surv_time, vital_status) ~ is_module_mutated, data = df.update)

					if (nchar(dframe) < 50) {
						title <- paste("Module:", dframe)
					} else {
						splits <- unlist(strsplit(dframe, split = ";"))
						if (length(splits) <= 20) {
							first <- paste(splits[1:10],  collapse = ";")
							second <- paste(splits[11:length(splits)], collapse = ";")
							title <- paste("Module: ", first, ";", "\n", second, sep = "")
						} else if (20 < length(splits) & length(splits) <= 30) {
							first <- paste(splits[1:10],  collapse = ";")
							second <- paste(splits[11:20], collapse = ";")
							third <- paste(splits[20:length(splits)], collapse = ";")
							title <- paste("Module: ", first, ";", "\n", second, ";", "\n", third, sep = "")
						} else {
						first <- paste(splits[1:10],  collapse = ";")
						second <- paste(splits[11:20], collapse = ";")
						third <- paste(splits[21:30], collapse = ";")
						fourth <- paste(splits[31:length(splits)], collapse = ";")
						title <- paste("Module: ", first, ";", "\n", second, ";", "\n", third, ";", "\n", fourth, sep = "")
						}
					}

					lst[[paste(cancer, dframe, sep = "")]] <- ggsurvplot(logrank.fit,
					   						data = df.update, 
					   						title = title,
					   						xlab = "Time (Days)",
					   						font.main = c(10, "bold"),
					   						font.caption = c(10),
											linetype = "strata", 
											conf.int = TRUE, 
											# pval = TRUE, pval.method = TRUE,
											risk.table = TRUE,
											legend = "bottom", legend.title = "Mutation Status", legend.labs = c("wildtype", "mutated"), 
											palette = c("#a1d0ff", "#ffc5a1"),
											caption = caption
					   			 )
		}
	}	

			res <- arrange_ggsurvplots(lst, print = FALSE, ncol = 1, nrow = 1)
			ggsave(survival.file, res)
			dev.off()
}

PatientsVSGenesMatrixPlot <- function(nonlinkers.added, mimp.modules, muts) {
	# Draw mutation matrices	

	# split mimp modules 
	mimp.modules.lst <- split(mimp.modules, f = mimp.modules$cancer)

	full.with.muts <- nonlinkers.added %>%
							mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
							group_by(cancer) %>% 
							arrange(coxph) %>% 
							unnest(List_of_patients) 

	mut.joined <- inner_join(full.with.muts, muts, by = c("List_of_patients" = "patient", "cancer" = "Study_Name")) 
	
	full.muts <- split(mut.joined, f = mut.joined$cancer) 

	lst = list()
	for (cancer in names(full.muts)) { 
		pvalue.arranged <- full.muts[[cancer]] %>% 
								select(Module, coxph) %>% 
								arrange(coxph) %>%
								unique()

		full.muts[[cancer]]$Module <- factor(full.muts[[cancer]]$Module, levels = pvalue.arranged$Module)
		full.muts.cancer <- split(full.muts[[cancer]], f = full.muts[[cancer]]$Module) 

		for (df in names(full.muts.cancer)) {
			full.muts.cancer[[df]]$Module <- as.character(full.muts.cancer[[df]]$Module)

			#variables
			intersect.genes <- unique(full.muts.cancer[[df]]$intersect)
			cancer.type <- unique(full.muts.cancer[[df]]$cancer)
			mod.genes <- unique(full.muts.cancer[[df]]$Module)
			patient.num <- length(unique(full.muts.cancer[[df]]$List_of_patients))
			coxph.value <- unique(full.muts.cancer[[df]]$coxph)
			mimp <- unique(full.muts.cancer[[df]]$mimp)

			print(paste(cancer.type, df))

			# molten <- full.muts.cancer[[df]] %>% 
			# 				ungroup() %>% 
			# 				select(List_of_patients, Module, gene) %>% 
			# 				mutate(mods = strsplit(Module, split = ";")) %>% 
			# 				unnest(mods) %>%
			# 				filter(gene %in% mod.genes) %>%
			# 				unique() %>% 
			# 				table() %>% 
			# 				melt()

			# molten$List_of_patients <- as.character(molten$List_of_patients)
			# molten$gene <- as.character(molten$gene)
			# molten[molten == 0] <- "wt"
			# molten[molten == 1] <- "mutated"

			redd <- full.muts.cancer[[df]] %>% 
				rename(gene = hugo_name) %>%
				select(Module, List_of_patients, gene, mimp, number_of_mimp_patients, cancer) %>% 
				mutate(Module = strsplit(Module, split = ";")) %>%
				unnest(Module) %>%  
				filter(gene == Module) %>% 
				unique()

			B <- dcast(List_of_patients ~ gene, data = redd, length)

			melted <- melt(table(redd))
			melted[melted == 0] <- "wt"
			melted[melted == 1] <- "mutated"
			melted %>% 
				rename(mutation_status = value) -> melted 
			melted$mutation_status <- factor(melted$mutation_status)
			melted[] <- lapply(melted, as.character)

			a <- unique(redd$gene)
			b <- unlist(strsplit(mod.genes, split = ";")) 
			if (length(setdiff(b, a)) != 0) {
				diff.genes <- setdiff(b, a) 
				C <- as.data.frame(replicate(length(diff.genes), data.frame(rep(0, nrow(B)))), stringsAsFactors=F)
				colnames(C) <- diff.genes
				tr <- cbind(B, C)
				tr <- melt(tr)
				tr$value <- as.character(tr$value)
				tr[tr == "0"] <- "wildtype"
				tr[tr == "1"] <- "mutated"
				tr$List_of_patients <- as.character(tr$List_of_patients)
				tr$variable <- as.character(tr$variable)
			} else {
				tr <- melt(B)
				tr[tr == "0"] <- "wildtype"
				tr[tr == "1"] <- "mutated"
				tr$List_of_patients <- as.character(tr$List_of_patients)
				tr$variable <- as.character(tr$variable)
				tr$value <- as.character(tr$value)
			}

			if (!cancer.type %in% names(mimp.modules.lst)) {
				joined <- tr %>% 
							mutate(mimp = FALSE)
			} else {
				if (mimp.modules.lst[[cancer.type]] %>%
					filter(Module == df) %>% nrow() == 0) {
				joined <- tr %>% 
								mutate(mimp = FALSE)
				} else {
					temp <- mimp.modules.lst[[cancer.type]] %>%
								filter(Module == df) %>% 
								mutate(patient = strsplit(patient, split = ","), patient_muts = strsplit(patient_muts, split = ";")) %>% 
								unnest(patient, patient_muts) %>% 
								mutate(gene = gsub(" .*", "", patient_muts), mimp = TRUE) %>%
								mutate(mimp = TRUE) %>% 
								select(gene, patient, mimp) %>%
								unique()

					joined <- tr %>% 
								left_join(temp, by = c("variable" = "gene", "List_of_patients" = "patient"))

					joined[is.na(joined)] <- FALSE
				}
			}

		names(sort(table(joined[,c(2,3)])[,1])) -> sorting 
		joined$variable <- factor(joined$variable, levels = sorting)

		names(sort(table(joined[,c(1,3)])[,1], decreasing = TRUE)) -> id.sorting
		joined$List_of_patients <- factor(joined$List_of_patients, levels = id.sorting)

		lst[[paste(cancer, df, sep = "_")]] <- ggplot(joined, aes(x = List_of_patients,  y = variable, fill = value)) + 
													geom_tile(colour = "white", size = 1.8) + 
													scale_y_discrete(expand =c(0,0)) + 
													scale_x_discrete(expand=c(0,0)) +
										   			geom_point(data = joined[joined$mimp, ], aes(size = factor(as.numeric(mimp), labels = c("mutated")))) +
													theme_grey(base_size = 18)+ 
													theme(plot.title = element_text(size = 40, face = "bold"),
														plot.caption = element_text(size = 35),
														axis.text = element_text(size = 30), 
														axis.title = element_text(size = 35),
														axis.ticks = element_line(size = 0.8),
														legend.title = element_text(size = 35),
														legend.text = element_text(size = 30),
														legend.key.size = unit(4, 'lines'),
														legend.position = "bottom",
														plot.background = element_blank(),
														panel.border = element_blank(), 
														axis.text.x = element_blank()) + 
													scale_fill_manual(values = c(wildtype = "#a1d0ff", mutated = "#ffc5a1"), guide = "legend") +
													coord_fixed() + 
													guides(fill = guide_legend(override.aes = list(shape = ''))) + 
										         	labs(x = "Patient",
														y = "Gene", 
														size = "mimp mutation status",
														fill = "gene mutation status",
														title = "Patient Mutation Table",
														caption = paste(cancer.type, intersect.genes, sep = "\n"))
			}
	}
	return(lst)
}

WriteCytoscapeNetworks <- function(nonlinkers.added, network, network.folder) {
	# Write out modules for cytoscape

	splitted <- split(nonlinkers.added, f = nonlinkers.added$cancer)
	arranged <- lapply(splitted, function(dframe) dframe %>% arrange(coxph))

	sorted <- bind_rows(arranged) # order 

	selected <- sorted %>%
					select(Module) %>% 
					mutate(mod.genes = strsplit(Module, split = ";")) %>% 
					unnest()

	joined <- inner_join(selected, network, by = c("mod.genes" = "V1"))

	# split by module 
	splitted.2 <- split(joined, f = joined$Module, drop = TRUE)

	# reorder
	ex <- splitted.2[sorted$Module]

	# counter for name
	i <- 1
	for (mod in names(ex)) {

		mod.dframe <- data.frame(genes = unlist(strsplit(mod, split = ";")), stringsAsFactors = FALSE)

		# remove non module genes in substrates 
		joined <- inner_join(ex[[mod]][, c(2,3)], mod.dframe, by = c("V2" = "genes")) %>% 
						filter(mod.genes != V2)

		# remove duplication
		new.dframe <- as.data.frame(unique(t(apply(joined, 1, sort))), stringsAsFactors = FALSE)

		write.table(new.dframe, file = paste(network.folder, "net", i, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

		i <- i + 1 
	}
}
############################################################################################################################################
mut.joined <- JoinModuleMutations(mimp.added.modules, muts)
nonlinkers.added <- AddNonLinker(mut.joined)
saveRDS(nonlinkers.added, modules.non.linker.added.file)

DrawSurvivalCurvesByCancerType(nonlinkers.added, muts, clins, survival.file) 

patient.matrices <- PatientsVSGenesMatrixPlot(nonlinkers.added, mimp.modules, muts) 
pdf(mutation.matrix.file, width = 70, height = 28, onefile = T)
		for (df in names(patient.matrices)) {
			print(patient.matrices[[df]])
	}

dev.off()

WriteCytoscapeNetworks(nonlinkers.added, network, network.folder)