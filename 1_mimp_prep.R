# mutation data prep for MIMP 

require(data.table); require(tidyr); require(dplyr); require(reshape2); 
require(hash); require(magrittr); require(ggpubr); require(ggsci); 
require(gProfileR); require(gridExtra); require(cowplot);
require(rvest); require(xml2)

populate_uniques <- function(unique_input, duplicated_input) {
	# requires a while loop to repeat function, probably too slow
	# return merged duplicated_input and unique_input by adding to entries in duplicated_input to comments section of unique_input 
    merged <- merge(unique_input, duplicated_input, by = c(colnames(duplicated_input)[-5]), all=TRUE)
    # only join together last two columns if not NA 
    merged$Tumor_Sample_Barcode.y[is.na(merged$Tumor_Sample_Barcode.y)] <- ""
    merged$Tumor_Sample_Barcode <- trimws(paste(merged$Tumor_Sample_Barcode.x, merged$Tumor_Sample_Barcode.y))
    merged[,c(5,6)] <- NULL
    merged

#     # until nrow(unique(duplicated_input)) == nrow(duplicated_input): 
#     popoulate_uniques(duplicated_input, )

#     helper <- function(duplicated_input) {

#     	# base case 

#     }
# }

populate_uniques1 <- function(duplicated_input) {
	# works but too slow 
	for (i in 1:nrow(duplicated_input)) {
		flag <- paste(duplicated_input[i, 1], duplicated_input[i, 2], duplicated_input[i, 3], duplicated_input[i, 4])
		index = which(paste(unique_input$Chromosome, unique_input$Start_Position, unique_input$Reference_Allele, unique_input$Tumor_Seq_Allele2) %in% flag)
		unique_input[index, 5] <- paste(unique_input[index, 5], duplicated_input[i, 5])
	}
	unique_input
}

helper <- function(row) {
	flag <- paste(row[1], row[2], row[3], row[4])
		index = which(paste(unique_input$Chromosome, unique_input$Start_Position, unique_input$Reference_Allele, unique_input$Tumor_Seq_Allele2) %in% flag)
		unique_input[index, 5] <- paste(unique_input[index, 5], row[5])
}


func <- function(duplicated_input, unique_input) {
	# works but too slow
	adply(duplicated_input, 1, function(row) {
	flag <- paste(row[1], row[2], row[3], row[4])
		index = which(paste(unique_input$Chromosome, unique_input$Start_Position, unique_input$Reference_Allele, unique_input$Tumor_Seq_Allele2) %in% flag)
		unique_input[index, 5] <- paste(unique_input[index, 5], row[5])
	})
	unique_input 
}

# populate_uniques1 <- function(duplicated_input, unique_input = unique_in) {
# 	apply(duplicated_input, 1, helper(duplicated_input))
# 	unique_input


# make it do something for each row so get rid of the for loop within the for loop 
# 	row.names(duplicated_input) <- NULL
# 	row.names(unique_input) <- NULL
# 	for (i in 1:nrow(duplicated_input)) {
# 		if (identical(duplicated_input[i, -5], unique_input[i, -5])) {
# 			unique_input[i, 5] = paste(unique_input[i, 5], duplicated_input[i, 5])
# 		}
# 	}
# 	unique_input


populate_uniques2 <- function(duplicated_input, unique_input) {
	# requires hash package 
	key <- do.call("paste", list(unique_input$Chromosome, unique_input$Start_Position, unique_input$Reference_Allele, unique_input$Tumor_Seq_Allele2))
	value <- c(unique_input$Tumor_Sample_Barcode)
	h <- hash(key, value)
	apply(duplicated_input, 1, function(row) {
		k <- paste(row[1], row[2], row[3], row[4])
		h[[k]] <- paste(h[[k]], row[5])
		})
	h
}

populate_uniques3 <- function(input) {
	# need to recombine hash tables because too big
	h <- hash()
	apply(input, 1, function(row) {
		k <- paste(row[1], row[2], row[3], row[4])
		if (!has.key(k, h)) {
			.set(h, k, row[5])
		} else {
			h[[k]] <- paste(h[[k]], row[5])
		} 
	})
	h
}

populate_uniques4 <- function(input) {
	# this one works
	aggregates <- aggregate(input, by = input[-6], paste)
	aggregates[, c(1:5, 11)]
}

wildtype_match <- function(mut_data, fasta_file) {
	nonmatched <- sapply(mut_data, function(row) {
		ref <- gsub(" [ABCDEFGHIJKLMNOPQRSTUVWXYZ][0123456789]+[ABCDEFGHIJKLMNOPQRSTUVWXYZ]", "", row)
		mut <- gsub("NM_[0123456789]+ ", "", row)
		if (length(fasta_file[which(names(fasta_file) %in% ref)]) == 0) {
			paste(ref, "NA")
		} else if (substr(fasta_file[which(names(fasta_file) %in% ref)], substr(mut, 2, nchar(mut)-1), substr(mut, 2, nchar(mut)-1)) != substr(mut, 1, 1)) {
			row
		} else if (substr(fasta_file[which(names(fasta_file) %in% ref)], substr(mut, 2, nchar(mut)-1), substr(mut, 2, nchar(mut)-1)) == substr(mut, 1, 1)) {
			"good"
		} else {
			"bad"
		}
	})
	which(!unlist(unname(nonmatched)) %in% "good")
}

refs <- sapply(dontmatches, function(row) {
		gsub(" [ABCDEFGHIJKLMNOPQRSTUVWXYZ][0123456789]+[ABCDEFGHIJKLMNOPQRSTUVWXYZ]", "", row) }
		)


###########################################################################################################################################
maf.file <- fread("/.mounts/labs/reimandlab/private/users/thuang/data/05-29-17/full_maf.bed") #3600963 lines
colnames(full_maf) = full_maf[1, ]
full_maf <- full_maf[-1, ] # full_maf.rds

# prepare annovar input 

tumour_id <- which(colnames(full_maf) == "Tumor_Sample_Barcode")
chromosome_num <- which(colnames(full_maf) == "Chromosome")
dna_start_pos <- which(colnames(full_maf) == "Start_Position")
dna_end_pos <- which(colnames(full_maf) == "End_Position")
wt_allele <- which(colnames(full_maf) == "Reference_Allele")
stopifnot(identical(full_maf$Reference_Allele, full_maf$Tumor_Seq_Allele1))
mut_allele <- which(colnames(full_maf) == "Tumor_Seq_Allele2")

annovar_input <- data.frame(full_maf[, c(chromosome_num, dna_start_pos, dna_end_pos, wt_allele, mut_allele, tumour_id)]) # raw_annovar_input.rds
annovar_input <- unique(annovar_input) # 3600963 -> 3451181
# unique_input <- annovar_input[which(duplicated(annovar_input[, -5]) == FALSE), ] # 3091454
# unique_input$Tumor_Sample_Barcode <- paste("comment:", unique_input$Tumor_Sample_Barcode, sep = " ") 

# # populate unique_input$Tumor_Sample_Barcode
# duplicated_input <- annovar_input[which(duplicated(annovar_input[, -5]) == TRUE),] # 359727
prepped_annovar_input <- populate_uniques4(annovar_input) # 3091453
# add back line with NA 
prepped_annovar_input <- rbind(prepped_annovar_input, data.frame("Chromosome" = 19, "Start_Position" = 7267392, "End_Position" = 7267392, "Reference_Allele" = "C", "Tumor_Seq_Allele2" = NA, "Tumor_Sample_Barcode" = "TCGA-13-0889-01A-01W-0420-08"))
# collapse tumour_id list to character string 
prepped_annovar_input[, 6] <- apply(prepped_annovar_input, 1, function(row) paste(unlist(row[6]), collapse = " ")) 
prepped_annovar_input[, 2] <- as.integer(prepped_annovar_input[, 2])
prepped_annovar_input[, 1] <- as.integer(prepped_annovar_input[, 1])
prepped_annovar_input[, 3] <- as.integer(prepped_annovar_input[, 3]) # final_annovar_input.rds
write.table(prepped_annovar_input, file = "prepped_annovar_input.avinput", sep = "\t", quote = F, col.names = F, row.names = F) # 3091454

# annovar
perl table_annovar.pl /.mounts/labs/reimandlab/private/users/thuang/annovar_test_cbw/annovar/prepped_annovar_input.avinput --outfile final_output /.mounts/labs/reimandlab/private/users/thuang/annovar_test_cbw/annovar/humandb/ -buildver hg19 -remove -otherinfo -protocol refGene -operation g -nastring .
# WARNING: A total of 356 sequences will be ignored due to lack of correct ORF annotation
# 31486 in final_output.refGene.invalid_input: 3       107491719       107491720       -       A       TCGA-EO-A3KX-01A-11D-A228-09
# 3091455 final_output.hg19_multianno.txt although 3091454 prepped_annovar_input.avinput

annovar_output <- read.delim("/.mounts/labs/reimandlab/private/users/thuang/data/06-05-17/final_output.hg19_multianno.txt", stringsAsFactors = F)
# take only missense ones from annovar 
unique(annovar_output$ExonicFunc.refGene)
 # [1] "."                          "frameshift deletion"       
 # [3] "stopgain"                   "unknown"                   
 # [5] "nonframeshift deletion"     "stoploss"                  
 # [7] "nonsynonymous SNV"          "synonymous SNV"            
 # [9] "frameshift substitution"    "nonframeshift substitution"

nonsynonymous_SNV <- which(annovar_output$ExonicFunc.refGene %in% "nonsynonymous SNV")
length(nonsynonymous_SNV) # 1580384 (1558984 from take 1)
missense <- annovar_output[nonsynonymous_SNV, ]
missense <- missense$AAChange.refGene #full_missense.rds
missense <- missense$AAChange.refGene

source("/.mounts/labs/reimandlab/private/users/thuang/bin/mimp_inscript2.R") 

muts <- "/.mounts/labs/reimandlab/private/users/thuang/data/06-05-17/longest_refs2.rds"
mut_data <- readRDS(muts) #1580384
# delete any lines with wrong structure
which(longest_refs2 %in% -Inf) # 88664  355083  665164 1186183
mut_data <- mut_data[-grep("NM_[0123456789]+ [ABCDEFGHIJKLMNOPQRSTUVWXYZ][0123456789]+[ABCDEFGHIJKLMNOPQRSTUVWXYZ]", mut_data, invert = TRUE)] #1580380
wildtype_match(mut_data, fa1) # make sure wt residues in mut_data match amino acid in fa1 , 856 ones that don't match or refseq not in fa1; mut_data_doesnt_match_fasta_wt_aa.rds
write.table(mut_data, file = "mut_data.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)