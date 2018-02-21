# do validation in pcawg dataset 

require(data.table); require(tidyr); require(dplyr); require(reshape2)
require(survival); require(survminer)
require(gProfileR); require(stringr)

# constants
coxph.threshold <- 0.15

# files 
pcawg.maf <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/filtered_pcawg.tsv")
pcawg.converts <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/pcawgConversion.tsv")
pcawg.surv <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/Jan26_PCAWG_clinical")

tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/modules_added_nonlinkers.rds")
tcga.clins <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_fixed_clin.rds")
tcga.muts <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_mutation_data.rds")
mimp.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-05-18/mimp_modules.rds")

# file names
scurves.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/pcawg_scurves.pdf"
tcga.pcawg.scurves.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/tcga_pacwg_scurves.pdf"
mutation.matrix.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/pcawg_mutation_matrices.pdf"
tcga.mutation.matrix.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-12-18/tcga_pcawg_mutation_matrices.pdf"
small.modules.filtered.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/pcawg_modules.rds"

FilterMAF <- function(pcawg.maf) {
	# Filter maf file for only missense" SNP

	filtered <- pcawg.maf %>% 
					filter(Variant_Classification == "Missense_Mutation" & Variant_Type == "SNP")
	return(filtered)
}

FilterTCGAOverlapPatients <- function(pcawg.filtered, pcawg.surv) {
	# filter overlapping tcga and pcawg patients

	joined <- pcawg.surv %>% 
					filter(nchar(submitted_donor_id.y) !=12) %>% # remove tcga patients
					select(icgc_donor_id, donor_vital_status, donor_survival_time, donor_interval_of_last_followup) %>%
					mutate(survival_time = pmax(donor_survival_time, donor_interval_of_last_followup, na.rm = T)) %>%
					filter(!is.na(donor_vital_status)) %>%
					filter(survival_time > 0) %>%
					filter(!is.na(survival_time)) %>%
					unique() %>% 
					inner_join(pcawg.filtered, by = c("icgc_donor_id" = "Donor_ID")) %>%
					select(-donor_interval_of_last_followup, -donor_survival_time)
	return(joined) 
}

MakePCAWGTCGACancerKey <- function(tcga.modules, pcawg.joined) {
	# make a key to convert cancer types between tcga and pcawg

	# make list with tcga cancer types as names
	lst <- rep("", length(unique(tcga.modules$cancer)))
	names(lst) <- unique(tcga.modules$cancer)
	
	# manual matching 
	lst[["Kidney_renal_clear_cell_carcinoma"]] <- "Kidney-RCC"
	lst[["Lymphoid_Neoplasm_Diffuse_Large_B-cell_Lymphoma"]] <- "Lymph-BNHL"
	lst[["Pancreatic_adenocarcinoma"]] <- "Panc-AdenoCA"
	lst[["Liver_hepatocellular_carcinoma"]] <- "Liver-HCC"

	return(lst)
}

GetPCAWGPatients <- function(tcga.modules, pcawg.joined, pcawg.tcga.key) {
	# return pcawg patients that have module mutations by cancer type
 	
 	# get dataframe of keys 
 	dframe <- as.data.frame(pcawg.tcga.key, stringsAsFactors=F)
 	dframe$tcga.cancer <- rownames(dframe)
 	row.names <- NULL
	isolates <- filter(dframe, pcawg.tcga.key != "")

	# join with pcawg keys 
 	joined <- inner_join(tcga.modules, isolates, by = c("cancer" = "tcga.cancer"))

 	# join with pcawg.joined 
 	joined.2 <- joined %>%
 					mutate(pcawg.tcga.key = strsplit(pcawg.tcga.key, split = ",")) %>%
 					unnest(pcawg.tcga.key) %>%
 					mutate(intersect.genes = strsplit(intersect, split = ";")) %>%
 					unnest(intersect.genes) %>%
 					inner_join(pcawg.joined, by = c("pcawg.tcga.key" = "Project_Code", "intersect.genes" = "Hugo_Symbol")) %>% # join together with pcawg patients
					select(Module, intersect, cancer, intersect.genes, pcawg.tcga.key, icgc_donor_id) %>%
					group_by(cancer, Module) %>%
					summarize(intersect.genes = paste(unique(intersect.genes), collapse = ";"), intersect = paste(unique(intersect), collapse = ";"), icgc_donor_id = paste(unique(icgc_donor_id), collapse = ",")) 

	return(joined.2)
}

SplitPCAWGModules <- function(pcawg.modules) {
	# split modules by cancer type

	splitted.modules <- split(pcawg.modules, pcawg.modules$cancer)

	return(splitted.modules)

}

SplitPCAWGPatients <- function(pcawg.joined, pcawg.tcga.key) {
	# format and split pcawg.surv 

	dframe <- as.data.frame(pcawg.tcga.key, stringsAsFactors=F)
 	dframe$tcga.cancer <- rownames(dframe)
 	row.names <- NULL
	isolates <- filter(dframe, pcawg.tcga.key != "")

	unnested.isolates <- isolates %>% 
								mutate(pcawg.tcga.key = strsplit(pcawg.tcga.key, split = ",")) %>% 
								unnest(pcawg.tcga.key)

	pcawg.patients <- pcawg.joined %>% 
							inner_join(unnested.isolates, by = c("Project_Code" = "pcawg.tcga.key")) %>%
							select(icgc_donor_id, donor_vital_status, survival_time, tcga.cancer) %>% 
							unique()

	splitted.pcawg.patients <- split(pcawg.patients, pcawg.patients$tcga.cancer)

	return(splitted.pcawg.patients)
}

h_GetPCAWGCoxphLogRank <- function(modules = splitted.modules, patients) {
	# Get coxph scores for modules in each cancer type 

	splitted <- split(modules, modules$Module)

	# format for coxph 
	split.rows <- list() 
	for (df in names(splitted)) {
		df_update <- splitted[[df]] %>%
						ungroup() %>%
						mutate(icgc_donor_id = strsplit(icgc_donor_id, split = ",")) %>% 
						select(icgc_donor_id, intersect.genes) %>%
						unnest(icgc_donor_id) %>%
						unique() %>%
						mutate(is_module_mutated = 1) %>%
						full_join(patients, by = "icgc_donor_id") %>%
						mutate(donor_vital_status = gsub("deceased", 1, donor_vital_status)) %>%
						mutate(donor_vital_status = gsub("alive", 0, donor_vital_status)) %>%
						mutate(donor_vital_status = as.numeric(donor_vital_status)) %>% 
						select(icgc_donor_id, is_module_mutated, donor_vital_status, survival_time)

					df_update$is_module_mutated[is.na(df_update$is_module_mutated)] <- 0
					split.rows[[df]] <- df_update
	}

	# calculate logrank 
	lst = list()
	for (df in names(split.rows)) {
		tryCatch(sdiff <- survdiff(Surv(survival_time, donor_vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "sdiff error")))
		sdiff <- survdiff(Surv(survival_time, donor_vital_status) ~ is_module_mutated, data = split.rows[[df]])
		p_logrank = 1-pchisq(sdiff$chisq, df=1)
		lst[[df]] <- p_logrank
	}

	chr <- as.numeric(as.character(lst))
	names(chr) <- names(lst)

	# calculate coxph
	lst_c <- list()
	lst_coeff <- list() # correlation
	for (df in names(split.rows)) {
		tryCatch(h0 <- coxph(Surv(survival_time, donor_vital_status) ~ 1, data = split.rows[[df]]), warning = function(w) print(paste(df, "h0 error")))
		h0 <- coxph(Surv(survival_time, donor_vital_status) ~ 1, data = split.rows[[df]])
		tryCatch(h1 <- coxph(Surv(survival_time, donor_vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "h1 error")))
		h1 <- coxph(Surv(survival_time, donor_vital_status) ~ is_module_mutated, data = split.rows[[df]])
		p_cox <- anova(h0, h1)[2,4] 
		lst_c[[df]] <- p_cox

		# correlation
		coeff <- h1[[1]]
		names(coeff) <- NULL
		lst_coeff[[df]] <- coeff 
	}

	chr_c <- as.numeric(as.character(lst_c))
	names(chr_c) <- names(lst_c)

	chr_c_corr <- as.numeric(as.character(lst_coeff))
	names(chr_c_corr) <- names(lst_coeff)

	dat <- data.frame(chr, chr_c, chr_c_corr) %>%
				tibble::rownames_to_column() %>%
				rename(module = rowname, logrank = chr, coxph = chr_c, coxph_corr = chr_c_corr) %>%
				inner_join(modules, by = c("module" = "Module"))
	return(dat)
}

RemoveNonsignificantCoxph <- function(test.results, coxph.threshold) {
	# Remove modules with coxph greater than 0.05

	bind_rows(test.results) %>% 
		filter(coxph <= coxph.threshold) -> filtered

return(filtered)		
}

FilterTestResults <- function(test.results) {
	# Filter out small modules by patients and by number of genes in modules

	small.modules.filtered <- bind_rows(test.results) %>% 
								rename(Module = module, calculated_logrank_pvalue = logrank, calculated_coxph_pvalue = coxph, List_of_patients = icgc_donor_id) %>%
								mutate(Number_patients = str_count(List_of_patients, ",") + 1) %>%
								filter(Number_patients >= 5) %>% # 19 left
								filter((str_count(intersect.genes, ";") + 1) >= 3) %>% 
								filter(!duplicated (intersect.genes)) # 18 left

	return(small.modules.filtered)
}

FormatClinResults <- function(pcawg.tcga.key, pcawg.joined) {
	# Format clinical data 

	# format pcawg.joined
	dframe <- as.data.frame(pcawg.tcga.key, stringsAsFactors=F)
	dframe$tcga.cancer <- rownames(dframe)
	row.names <- NULL
	isolates <- dframe %>% 
					filter(pcawg.tcga.key != "") %>% 
					mutate(pcawg.tcga.key = strsplit(pcawg.tcga.key, split = ",")) %>% 
					unnest(pcawg.tcga.key)

	clin <- unique(inner_join(pcawg.joined, isolates, by = c("Project_Code" = "pcawg.tcga.key")))

	return(clin)
}

FilteredSurvivalCurvesByCancerType <- function(small.modules.filtered, pcawg.joined, pcawg.tcga.key, file.name = scurves.file) {
	# Filter out small modules and draw survival curves 

	filtered.muts <- split(small.modules.filtered, f = small.modules.filtered$cancer)

	clin.lst <- split(clin, f = clin$tcga.cancer)

	lst <- list()
	for (cancer in names(filtered.muts)) { 
		filtered.muts[[cancer]] %>% 
			select(intersect.genes, calculated_coxph_pvalue) %>% 
			arrange(calculated_coxph_pvalue) %>% 
			unique() -> tt 

		filtered.muts[[cancer]]$intersect.genes <- factor(filtered.muts[[cancer]]$intersect.genes, levels = tt$intersect.genes)
		filtered.muts.cancer <- split(filtered.muts[[cancer]], f = filtered.muts[[cancer]]$intersect.genes) 

		for (df in names(filtered.muts.cancer)) {
			filtered.muts.cancer[[df]]$intersect.genes <- as.character(filtered.muts.cancer[[df]]$intersect.genes)
			intersect.genes <- unique(filtered.muts.cancer[[df]]$intersect.genes)
			coxph.value <- formatC(unique(filtered.muts.cancer[[df]]$calculated_coxph_pvalue), format = "e", digits = 2)
			logrank.value <- formatC(unique(filtered.muts.cancer[[df]]$calculated_logrank_pvalue), format = "e", digits = 2)


			filtered.muts.cancer[[df]] %>%
				mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
				select(intersect.genes, List_of_patients) %>%
				unnest(List_of_patients) %>%
				unique() %>%
				mutate(is_module_mutated = 1) %>%
				full_join(clin.lst[[cancer]], by = c("List_of_patients" = "icgc_donor_id")) %>%
				mutate(donor_vital_status = gsub("deceased", 1, donor_vital_status)) %>%
				mutate(donor_vital_status = gsub("alive", 0, donor_vital_status)) %>%
				mutate(donor_vital_status = as.numeric(donor_vital_status)) %>% 
				select(List_of_patients, is_module_mutated, donor_vital_status, survival_time) %>% 
				unique() -> df_update

				df_update$is_module_mutated[is.na(df_update$is_module_mutated)] <- 0

				caption <- paste(paste("logrank p = ",logrank.value, "; coxph p = ", coxph.value, sep = ""), cancer, sep = "\n")

				logrank.fit <- survfit(Surv(survival_time, donor_vital_status) ~ is_module_mutated, data = df_update)

					if (nchar(df) < 50) {
						title <- paste("Module:", df)
					} else {
						splits <- unlist(strsplit(df, split = ";"))
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

					lst[[paste(cancer, df, sep = "")]] <- ggsurvplot(logrank.fit,
					   						data = df_update, 
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
	ggsave(file.name, res)
}	

JoinPatientMutations <- function(small.modules.filtered, clin) {
	# Join patient mutations

	joined <- small.modules.filtered %>% 
					mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>%
					unnest(List_of_patients) %>% 
					inner_join(clin, by = c("List_of_patients" = "icgc_donor_id")) %>%
					unique()

	joined.lst <- split(joined, joined$cancer)
	
	return(joined.lst)
}

PatientsVSGenesMatrixPlot <- function(full.muts = patient.mutations.lst) {

	lst <- list()
	for (cancer in names(full.muts)) { 
		full.muts[[cancer]] %>% 
			select(intersect.genes, calculated_coxph_pvalue) %>% 
			arrange(calculated_coxph_pvalue) %>% 
			unique() -> tt 

		full.muts[[cancer]]$intersect.genes <- factor(full.muts[[cancer]]$intersect.genes, levels = tt$intersect.genes)
		full.muts.cancer <- split(full.muts[[cancer]], f = full.muts[[cancer]]$intersect.genes) 

		for (df in names(full.muts.cancer)) {
			full.muts.cancer[[df]]$intersect.genes <- as.character(full.muts.cancer[[df]]$intersect.genes)

			#variables
			cancer.type <- unique(full.muts.cancer[[df]]$cancer)
			mod.genes <- unique(full.muts.cancer[[df]]$intersect.genes)
			patient.num <- length(unique(full.muts.cancer[[df]]$List_of_patients))
			coxph.value <- unique(full.muts.cancer[[df]]$calculated_coxph_pvalue)

			print(paste(cancer.type, df))


			full.muts.cancer[[df]] %>% 
				select(intersect.genes, List_of_patients, Hugo_Symbol, cancer) %>% 
				mutate(intersect.genes = strsplit(intersect.genes, split = ";")) %>%
				unnest(intersect.genes) %>%  
				filter(Hugo_Symbol == intersect.genes) %>% 
				unique() -> redd

			B <- dcast(List_of_patients ~ Hugo_Symbol, data = redd, length)

			melted <- melt(table(redd))
			melted[melted == 0] <- "wt"
			melted[melted == 1] <- "mutated"
			melted %>% 
				rename(mutation_status = value) -> melted 
			melted$mutation_status <- factor(melted$mutation_status)
			melted[] <- lapply(melted, as.character)

			a <- unique(redd$Hugo_Symbol)
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

				tr %>% 
					mutate(mimp = FALSE) -> joined

			names(sort(table(joined[,c(2,3)])[,1])) -> sorting 
			joined$variable <- factor(joined$variable, levels = sorting)

			names(sort(table(joined[,c(1,3)])[,1], decreasing = TRUE)) -> id.sorting
			joined$List_of_patients <- factor(joined$List_of_patients, levels = id.sorting)

		lst[[paste(cancer, df, sep = "~")]] <- ggplot(joined, aes(x = List_of_patients,  y = variable, fill = value)) + 
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
					caption = paste(cancer.type))
		}
	}
	return(lst)
}

DrawTCGAMatchingSurvivalCurves <- function(small.modules.filtered, tcga.modules, tcga.clins, file.name = tcga.pcawg.scurves.file) {
	# Filter out small modules and draw survival curves 

	selected.small.modules <- small.modules.filtered %>% 
									select(Module, calculated_coxph_pvalue, cancer, intersect.genes, coxph_corr) %>% 
									rename(pcawg_coxph = calculated_coxph_pvalue, pcawg_coxph_corr = coxph_corr)

	intersected.modules <- inner_join(selected.small.modules, tcga.modules, by = c("cancer", "Module"))

	filtered.muts <- split(intersected.modules, f = intersected.modules$cancer)

	clin.lst <- split(tcga.clins, f = tcga.clins$Study_Name)

	lst <- list()
	for (cancer in names(filtered.muts)) { 
		filtered.muts[[cancer]] %>% 
			arrange(pcawg_coxph) %>% 
			unique() -> tt 

		# make the orders static
		filtered.muts[[cancer]]$Module <- factor(filtered.muts[[cancer]]$Module, levels = tt$Module)
		filtered.muts.cancer <- split(filtered.muts[[cancer]], f = filtered.muts[[cancer]]$Module) 

		for (df in names(filtered.muts.cancer)) {
			filtered.muts.cancer[[df]]$Module <- as.character(filtered.muts.cancer[[df]]$Module)
			Module <- unique(filtered.muts.cancer[[df]]$Module)
			intersect.genes <- unique(filtered.muts.cancer[[df]]$intersect.genes)
			coxph.value <- formatC(unique(filtered.muts.cancer[[df]]$coxph), format = "e", digits = 2)
			logrank.value <- formatC(unique(filtered.muts.cancer[[df]]$logrank), format = "e", digits = 2)


			df_update <- filtered.muts.cancer[[df]] %>%
								mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
								select(Module, List_of_patients) %>%
								unnest(List_of_patients) %>%
								unique() %>%
								mutate(is_module_mutated = 1) %>%
								full_join(clin.lst[[cancer]], by = c("List_of_patients" = "patient")) %>%
								mutate(vital_status = gsub("DECEASED", 1, vital_status)) %>%
								mutate(vital_status = gsub("ALIVE", 0, vital_status)) %>%
								mutate(vital_status = as.numeric(vital_status)) %>% 
								select(List_of_patients, is_module_mutated, vital_status, surv_time) %>% 
								unique() 

				df_update$is_module_mutated[is.na(df_update$is_module_mutated)] <- 0

				caption <- paste(paste("logrank p = ",logrank.value, "; coxph p = ", coxph.value, sep = ""), intersect.genes, cancer, sep = "\n")

				logrank.fit <- survfit(Surv(surv_time, vital_status) ~ is_module_mutated, data = df_update)

					if (nchar(df) < 50) {
						title <- paste("Module:", df)
					} else {
						splits <- unlist(strsplit(df, split = ";"))
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

					lst[[paste(cancer, df, sep = "")]] <- ggsurvplot(logrank.fit,
					   						data = df_update, 
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
	ggsave(file.name, res)
}	


DrawTCGAPatientsVSGenesMatrixPlot <- function(small.modules.filtered, tcga.modules, tcga.muts) {
	# Draw mutation matrices	

	# split mimp modules 
	selected.small.modules <- small.modules.filtered %>% 
									select(Module, calculated_coxph_pvalue, cancer, intersect.genes, coxph_corr) %>% 
									rename(pcawg_coxph = calculated_coxph_pvalue, pcawg_coxph_corr = coxph_corr)

	intersected.modules <- inner_join(selected.small.modules, tcga.modules, by = c("cancer", "Module"))

	# tcga.modules.lst <- split(intersected.modules, f = intersected.modules$cancer)

	# tcga.modules.lst <- split(tcga.modules, f = tcga.modules$cancer)

	full.with.muts <- intersected.modules %>%
							mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
							group_by(cancer) %>% 
							arrange(pcawg_coxph) %>% 
							unnest(List_of_patients) 

	mut.joined <- inner_join(full.with.muts, tcga.muts, by = c("List_of_patients" = "patient", "cancer" = "Study_Name")) 
	
	full.muts <- split(mut.joined, f = mut.joined$cancer) 
	mimp.modules.lst <- split(mimp.modules, f = mimp.modules$cancer)

	lst = list()
	for (cancer in names(full.muts)) { 
		pvalue.arranged <- full.muts[[cancer]] %>% 
								ungroup() %>%
								select(Module, pcawg_coxph) %>%
								arrange(pcawg_coxph) %>%
								unique()

		full.muts[[cancer]]$Module <- factor(full.muts[[cancer]]$Module, levels = pvalue.arranged$Module)
		full.muts.cancer <- split(full.muts[[cancer]], f = full.muts[[cancer]]$Module) 

		for (df in names(full.muts.cancer)) {
			full.muts.cancer[[df]]$Module <- as.character(full.muts.cancer[[df]]$Module)

			#variables
			intersect.genes <- unique(full.muts.cancer[[df]]$intersect.genes)
			cancer.type <- unique(full.muts.cancer[[df]]$cancer)
			mod.genes <- unique(full.muts.cancer[[df]]$Module)
			patient.num <- length(unique(full.muts.cancer[[df]]$List_of_patients))
			coxph.value <- unique(full.muts.cancer[[df]]$coxph)
			mimp <- unique(full.muts.cancer[[df]]$mimp)

			print(paste(cancer.type, df))

			redd <- full.muts.cancer[[df]] %>% 
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

################################################################################################################
pcawg.filtered <- FilterMAF(pcawg.maf)
pcawg.joined <- FilterTCGAOverlapPatients(pcawg.filtered, pcawg.surv) #414 unique donor ids
pcawg.tcga.key <- MakePCAWGTCGACancerKey(tcga.modules, pcawg.joined)
pcawg.modules <- GetPCAWGPatients(tcga.modules, pcawg.joined, pcawg.tcga.key)
splitted.modules <- SplitPCAWGModules(pcawg.modules) 
splitted.patients <- SplitPCAWGPatients(pcawg.joined, pcawg.tcga.key)

test.results <- list() # 28
	for (i in names(splitted.modules)) {
		result <- h_GetPCAWGCoxphLogRank(splitted.modules[[i]], splitted.patients[[i]])
		test.results[[i]] <- result
	} 

small.modules.filtered <- FilterTestResults(test.results)
saveRDS(small.modules.filtered, small.modules.filtered.file)

clin <- FormatClinResults(pcawg.tcga.key, pcawg.joined)

# pcawg survival curves
FilteredSurvivalCurvesByCancerType(small.modules.filtered, pcawg.joined, pcawg.tcga.key)
dev.off()

coxphfiltered.results <- RemoveNonsignificantCoxph(test.results, coxph.threshold) # 3 remaining
patient.mutations.lst <- JoinPatientMutations(small.modules.filtered, clin)

# pcawg mutation matrices
patient.matrices <- PatientsVSGenesMatrixPlot()

pdf(mutation.matrix.file, width = 70, height = 28, onefile = T)
		for (df in names(patient.matrices)) {
			print(patient.matrices[[df]])
	}

dev.off()

# matching tcga survival curves 
DrawTCGAMatchingSurvivalCurves(small.modules.filtered, tcga.modules, tcga.clins, file.name = tcga.pcawg.scurves.file)
dev.off()

# matching tcga mutation matrices
tcga.patient.matrices <- DrawTCGAPatientsVSGenesMatrixPlot(small.modules.filtered, tcga.modules, tcga.muts)

pdf(tcga.mutation.matrix.file, width = 70, height = 28, onefile = T)
		for (df in names(tcga.patient.matrices)) {
			print(tcga.patient.matrices[[df]])
	}

dev.off()

