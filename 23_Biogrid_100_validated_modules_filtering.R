# filter hypermodules results

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(survival); require(survminer)

options(warn = 1)

# constants
hypermodules.dir <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/hypermodules_biogrid_100/"
hypermodules.prefix <- "hypermodules_100_" # number of permutations
patient.threshrold <- 5 # threshold for number of patients in a module
coxph.threshold <- 0.05

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/"
dir.create(loc)

# files
clins <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_fixed_clin.rds")

# file names
nosubset.results.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/03-05-18/nosubset_hypermodules_results.rds"

ReadHypermodulesResults <- function(hypermodules.dir, hypermodules.prefix) {
	# read in hypermodules results and remove cancers without results
	setwd(hypermodules.dir)
	file.paths <- list.files()

	lst <- lapply(file.paths, function(fil) {
		if (!file.size(fil) == 0) {
			fread(fil, skip = 3, verbose = TRUE)
			}
		})

	nam <- gsub(hypermodules.prefix, "", file.paths)
	names(lst) <- (gsub(".txt", "", nam)) 

	# remove NULLs
	lst[sapply(lst, is.null)] <- NULL

	for (i in names(lst)) {
		if (ncol(lst[[i]]) == 2) {
			lst[[i]] <- NULL
		}
	}

	binded <- bind_rows(lst, .id = "cancer")

	return(binded)
}

RemoveSmallModules <- function(hypermodules.results, patient.threshold) {
	# remove modules that have less than patient.threshold patients
	results <- hypermodules.results %>%
					filter(Number_patients >= 5) 
	return(results)
}

SplitHypermodulesResults <- function(hypermodules.results) {
	results.lst <- split(hypermodules.results, hypermodules.results$cancer)
	return(results.lst)
}

SplitClinicalResults <- function(clins, cancer.names = names(hypermodules.results.lst)) {
	clins.ab <- clins %>% 
					filter(Study_Name %in% cancer.names)

	splitted <- split(clins.ab, clins.ab$Study_Name)
	return(splitted)
}

h_CalculateLogRankCoxph <- function(result, clin) {
	# returns logrank and coxph for each module in a cancer type 
	splitted <- split(result, result$Module)
	split.rows <- list()
	for (df in names(splitted)) {
		df_update <- splitted[[df]] %>%
						mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
						select(Seed, Module, List_of_patients) %>%
						unnest(List_of_patients) %>%
						unique() %>%
						mutate(is_module_mutated = 1) %>%
						full_join(clin, by = c("List_of_patients" = "patient")) %>%
						mutate(vital_status = gsub("DECEASED", 1, vital_status)) %>%
						mutate(vital_status = gsub("ALIVE", 0, vital_status)) %>%
						mutate(vital_status = as.numeric(vital_status)) %>% 
						select(List_of_patients, is_module_mutated, vital_status, surv_time)

					df_update$is_module_mutated[is.na(df_update$is_module_mutated)] <- 0
					split.rows[[df]] <- df_update
	}
	# calculate logrank 
	lst = list()
	for (df in names(split.rows)) {
		tryCatch(sdiff <- survdiff(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "sdiff error")))
		sdiff <- survdiff(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]])
		p_logrank = 1-pchisq(sdiff$chisq, df=1)
		lst[[df]] <- p_logrank
	}

	chr <- as.numeric(as.character(lst))
	names(chr) <- names(lst)

	# calculate coxph & survival correlation
	lst_c <- list()
	lst_coeff <- list() # correlation
	for (df in names(split.rows)) {
		tryCatch(h0 <- coxph(Surv(surv_time, vital_status) ~ 1, data = split.rows[[df]]), warning = function(w) print(paste(df, "h0 error")))
		h0 <- coxph(Surv(surv_time, vital_status) ~ 1, data = split.rows[[df]])
		tryCatch(h1 <- coxph(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "h1 error")))
		h1 <- coxph(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]])
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
				inner_join(result, by = c("module" = "Module"))
	return(dat)
}


RemoveNonsignificantCoxph <- function(test.results, coxph.threshold) {
	# Remove modules with coxph greater than 0.05

	filtered <- bind_rows(test.results) %>% 
					filter(coxph <= coxph.threshold) 

return(filtered)		
}

RemoveRedundantPatientModules <- function(coxphfiltered.results) {
	# Remove modules that are subsets of larger modules

	names(coxphfiltered.results)[1] <- "Module" # rename
	coxphfiltered.results <- split(coxphfiltered.results, f = coxphfiltered.results$cancer)

	filtered.results <- list()

	for (result in names(coxphfiltered.results)) {
		splitted <- coxphfiltered.results[[result]] %>% 
						mutate(List_of_patients = strsplit(List_of_patients, split = ","))

		patient.list <- splitted$List_of_patients
		names(patient.list) <- coxphfiltered.results[[result]]$Module

		sorted <- character()
		for (lst in names(patient.list)) {
			ss <- sort(patient.list[[lst]])
			sorted[[lst]] <- paste(ss, collapse = ",")
		}

		splitted %>% 
			mutate(List_of_patients = sorted) %>% 
			mutate(List_of_patients = strsplit(List_of_patients, split = ",")) -> sorted.splitted

		sorted.test <- sorted.splitted$List_of_patients
		names(sorted.test) <- sorted.splitted$Module

		uniques <- unique(sorted.test)
			if (compare(sorted.test, uniques, allowAll = TRUE)[1] == FALSE) {
					dups <- sorted.test[names(sorted.test[which(sorted.test %in% sorted.test[duplicated(sorted.test)])])]
				} else {
					dups <- list()
				}
			
		sortings <- sorted.test
		for (entry in names(sorted.test)) {
			 x <- sorted.test[[entry]]
			 deleted <- sortings
			 deleted[[entry]] <- NULL
			 for (y in names(deleted)) {
			 	if (all(deleted[[y]] %in% x)) {
			 		sortings[[y]] <- NULL
			 	}
			}
		}

		new.lst <- c(sortings, dups)

		screened.list <- list()
		for (mod in names(new.lst)) {
		screened.list[[mod]] <- paste(new.lst[[mod]], collapse = ",")
		}

		new.dframe <- data.frame(Module = as.character(names(screened.list)), stringsAsFactors = F) %>% 
						inner_join(coxphfiltered.results[[result]], by = "Module")

		new.dframe <- new.dframe[, c("Module", "logrank", "coxph", "cancer", "coxph_corr", "Number_patients", "List_of_patients")]

		filtered.results[[result]] <- new.dframe
	}
	binded.results <- bind_rows(filtered.results)

	return(binded.results)	
}

GetSurvivalCorrelations <- function(nosubset.results) {
	# Get number of postiive and negative survival correlations

	wt.dobetter <- sum(nosubset.results$coxph_corr > 0)
	mut.dobetter <- sum(nosubset.results$coxph_corr < 0)

	nums <- list(wt.dobetter, mut.dobetter)
	names(nums) <- c("wt_dobetter", "mut_dobetter")

	return(nums)
}


########################################################################################################
hypermodules.results <- ReadHypermodulesResults(hypermodules.dir, hypermodules.prefix)
removed.results <- RemoveSmallModules(hypermodules.results, patient.threshold)

hypermodules.results.lst <- SplitHypermodulesResults(removed.results)
clins.lst <- SplitClinicalResults(clins)

test.results <- list()
	for (i in names(hypermodules.results.lst)) {
		result <- h_CalculateLogRankCoxph(hypermodules.results.lst[[i]], clins.lst[[i]])
		test.results[[i]] <- result
	}

coxphfiltered.results <- RemoveNonsignificantCoxph(test.results, coxph.threshold)
nosubset.results <- RemoveRedundantPatientModules(coxphfiltered.results)
saveRDS(nosubset.results, nosubset.results.file)

survival.correlations <- GetSurvivalCorrelations(nosubset.results)



